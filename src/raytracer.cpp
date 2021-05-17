// The JSON library allows you to reference JSON arrays like C++ vectors and JSON sceneObjects like C++ maps.

#include "raytracer.h"

#include <iostream>
#include <fstream>
#include <string>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/transform.hpp>

#include "json.hpp"

using json = nlohmann::json;

const char *PATH = "scenes/";

double fov = 60;
colour3 background_colour(0, 0, 0);

json scene;
json sceneLights;
json sceneObjects;

int numObjects = 0;
int numLights = 0;
const int maxObjects = 25;
const int maxTriangles = 500;
struct Object* objects[maxObjects];

struct Light* lights[maxObjects];
std::vector<float>* sceneAmbient;
const float tolerance = 0.0001;
const int maxReflect = 8;
const float AIR = 1.0002926;
const int TriGroupSize = 4;
const int areaSamples = 64;
const int glossSamples = 16;

struct BVHnode 
{
	int mesh = -1; //parent mesh, to get material properties
	bool isLeaf = false;

	int numTris = 0;//in this bounding volume

	point3 bvMin = point3(0, 0, 0); //bounding volume coords
	point3 bvMax = point3(0, 0, 0);

	int rangeMin = 0;//which triangles in the mesh
	int rangeMax = 0;

	BVHnode* left;
	BVHnode* right;

};

struct Object 
{
	bool isSphere = false;
	bool isPlane = false;
	bool isMesh = false;

	point3 position = point3(0, 0, 0);
	float radius = 0;
	point3 normal = point3(0, 0, 0);

	//material properties
	bool hasAmbient = false;
	bool hasDiffuse = false;
	bool hasAlbedo = false;
	bool hasRoughness = false;
	bool hasSpecular = false;
	bool isReflective = false;
	bool isTransmissive = false;
	bool isRefractive = false;
	bool isGlossy = false;

	point3 ambient = point3(0, 0, 0);
	point3 diffuse = point3(0, 0, 0);
	point3 specular = point3(0, 0, 0);
	point3 transmissive = point3(0, 0, 0);
	point3 reflective = point3(0,0,0);

	float shininess = 0;
	float refraction = 0;
	float gloss = 0;
	float albedo = 0;
	float roughness = 0;
	int numTriangles = 0;

	point3 triangles[maxTriangles][3] = {};

	
	point3 centerMass = point3(0, 0, 0);

	BVHnode* bvh;
};

struct Light
{
	bool isPoint = false;
	bool isDirectional = false;;
	bool isSpot = false;
	bool isAmbient = false;
	bool isArea = false;
	bool hasPosition = false;

	point3 position = point3(0, 0, 0);
	point3 color = point3(0, 0, 0);
	point3 direction = point3(0, 0, 0);

	float radius = 0;
	float cutoff = 0;
};

float nearT(point3 e, point3 d, int &obj, point3 &normal,bool pick);
point3 calcLighting(point3 d, point3 e, point3 p, point3 N, int material, int level, bool pick);
bool traverseBoundingVolumes(point3 e, point3 d, BVHnode* bvh, float &tmin, int& tri, bool pick);

/****************************************************************************/

// here are some potentially useful utility functions



json find(json &j, const std::string key, const std::string value) {
	json::iterator it;
	for (it = j.begin(); it != j.end(); ++it) {
		if (it->find(key) != it->end()) {
			if ((*it)[key] == value) {
				return *it;
			}
		}
	}
	return json();
}

glm::vec3 vector_to_vec3(const std::vector<float> &v) {
	return glm::vec3(v[0], v[1], v[2]);
}

void findOrthonormalVectors(point3 w, point3 &u, point3 &v)
{

	point3 x = point3(-w[2],w[0],w[1]);

	u = glm::cross(w, x);
	v = glm::cross(w, u);



	u = glm::normalize(u);
	v = glm::normalize(v);

	/*point3 k = point3(0, 0, 0);
	//orthonormal basis
	if (w.x > w.y && w.x > w.z) { k = point3(1, 0, 0); }
	else if (w.y > w.x && w.y > w.z) { k = point3(0, 1, 0); }
	else { k = point3(0, 0, 1); }

	u = glm::normalize((glm::cross(w, k)));
	v = glm::normalize((glm::cross(w, u)));*/
}

void transformMesh(int mesh)
{
	
	glm::vec4 center = glm::vec4(0, 0, 0, 1);
	center[0] = objects[mesh]->centerMass[0];
	center[1] = objects[mesh]->centerMass[1];
	center[2] = objects[mesh]->centerMass[2];

	point3 rotationAxis = point3(0, 1, 0.5);
	float angle = glm::radians(90.0f);
	glm::mat4x4 rotationMatrix = glm::rotate(angle, rotationAxis);

	 
	//translate to origin to rotate and scale properly
	glm::mat4x4 Otranslate =
		glm::translate(
			glm::mat4(),
			glm::vec3(center.x, center.y, center.z));
	glm::mat4x4 OinverseTranslate = glm::inverse(Otranslate); // to bring it back to where it was

	//translate
	glm::mat4x4 translate =
		glm::translate(
			glm::mat4(),
			glm::vec3(1, 0.0, 0.0));
	glm::mat4x4 inverseTranslate = glm::inverse(translate);

	//scale
	glm::mat4x4 scale =
		glm::scale(
			glm::mat4(),
			glm::vec3(0.5, 0.7, 0.5));


	//the complete transformation for each point
	glm::mat4x4 transform = translate * Otranslate *rotationMatrix * scale * OinverseTranslate;

	for (int i = 0; i < objects[mesh]->numTriangles;i++) 
	{
		for (int j = 0; j < 3; j++)
		{

			//convert to vec4
			glm::vec4 point = glm::vec4(0, 0, 0, 1);
			point[0] = objects[mesh]->triangles[i][j][0];
			point[1] = objects[mesh]->triangles[i][j][1];
			point[2] = objects[mesh]->triangles[i][j][2];

			point = transform * point;

			objects[mesh]->triangles[i][j][0] = point[0];
			objects[mesh]->triangles[i][j][1] = point[1];
			objects[mesh]->triangles[i][j][2] = point[2];		
		}
	}
}

struct BVHnode* createBoundingVolume(int mesh, int start, int end)
{

	float minx = 100, miny = 100;//, minz = 100;
	float maxx = -100, maxy = -100;//, maxz = -100;

	float maxz = -100;
	float minz = 100;

	BVHnode* bvh = new BVHnode;

	bvh->mesh = mesh;
	bvh->rangeMax = end;
	bvh->rangeMin = start;
	bvh->numTris = end - start + 1;
	//std::cout << "\nRangeMin: " << bvh->rangeMin << "\tRangeMax: " << bvh->rangeMax;
	

	if(end == start)
	{
		bvh->isLeaf = true;
		//std::cout << "\nCREATED A LEAF NODE";

		for (int j = 0; j < 3;j++)
		{
			float x = objects[mesh]->triangles[start][j][0];
			float y = objects[mesh]->triangles[start][j][1];
			float z = objects[mesh]->triangles[start][j][2];

			maxx = std::max(x, maxx);
			maxy = std::max(y, maxy);
			maxz = std::max(z, maxz); // since - z is farther

			minx = std::min(x, minx);
			miny = std::min(y, miny);
			minz = std::min(z, minz); // since + z is closer

			bvh->bvMin = point3(minx, miny, minz);
			bvh->bvMax = point3(maxx, maxy, maxz);

			
		}
		/*std::cout << "\nMesh Minx: " << minx;
		std::cout << "\tMesh Maxx: " << maxx;
		std::cout << "\nMesh Miny: " << miny;
		std::cout << "\tMesh Maxy: " << maxy;
		std::cout << "\nMesh MinZ: " << minz;
		std::cout << "\tMesh MaxZ: " << maxz;
		std::cout << "\nNum Tris in this BV: " << bvh->numTris;*/
	}
	else
	{
		for (int i = start; i <= end;i++)
		{
			for (int j = 0; j < 3;j++)
			{
				float x = objects[mesh]->triangles[i][j][0];
				float y = objects[mesh]->triangles[i][j][1];
				float z = objects[mesh]->triangles[i][j][2];

				maxx = std::max(x, maxx);
				maxy = std::max(y, maxy);
				maxz = std::max(z, maxz); 

				minx = std::min(x, minx);
				miny = std::min(y, miny);
				minz = std::min(z, minz); 
			}
		}

		

		bvh->bvMin = point3(minx, miny, minz);
		bvh->bvMax = point3(maxx, maxy, maxz);

		/*std::cout << "\nMesh Minx: " << minx;
		std::cout << "\tMesh Maxx: " << maxx;
		std::cout << "\nMesh Miny: " << miny;
		std::cout << "\tMesh Maxy: " << maxy;
		std::cout << "\nMesh MinZ: " << minz;
		std::cout << "\tMesh MaxZ: " << maxz;
		std::cout << "\nNum Tris in this BV: " << bvh->numTris;*/

		if (!bvh->isLeaf)
		{
			bvh->left = createBoundingVolume(mesh, start, (start + end) / 2);
			bvh->right = createBoundingVolume(mesh, 1+((start + end) / 2), end);
		}
	}
	return bvh;
}

void calcCenterMass(struct Object* mesh)
{
	float minx = 100, miny = 100;
	float maxx = -100, maxy = -100;

	float maxz = -100;
	float minz = 100;

	point3 totals = point3(0, 0, 0);

	for (int i = 0; i < mesh->numTriangles;i++)
	{
		point3 triavg = point3(0, 0, 0);
		for (int j = 0; j < 3;j++)
		{
			float x = mesh->triangles[i][j][0];
			float y = mesh->triangles[i][j][1];
			float z = mesh->triangles[i][j][2];

			maxx = std::max(x, maxx);
			maxy = std::max(y, maxy);
			maxz = std::max(z, maxz); 

			minx = std::min(x, minx);
			miny = std::min(y, miny);
			minz = std::min(z, minz);
			
			triavg[0] += x;
			triavg[1] += y;
			triavg[2] += z;

		}

		totals[0] += triavg[0] / 3;
		totals[1] += triavg[1] / 3;
		totals[2] += triavg[2] / 3;
	}

	mesh->centerMass = point3(totals.x / mesh->numTriangles, totals.y / mesh->numTriangles, totals.z / mesh->numTriangles);
	//std::cout << "\ncenterMass: (" << mesh->centerMass[0] << "," << mesh->centerMass[1] << "," << mesh->centerMass[2] << ")";
}



struct Object* createObject(json obj)
{
	 Object* object = new Object;

	//set type bools
	if (obj["type"] == "sphere") { object->isSphere = true; object->isPlane = false; object->isMesh = false; }
	else if (obj["type"] == "plane") { object->isPlane = true; object->isSphere = false; object->isMesh = false; }
	else if (obj["type"] == "mesh") { object->isMesh = true; object->isSphere = false; object->isPlane = false; }
	
	if (obj.find("radius") != obj.end()) { object->radius = (float)obj["radius"]; }
	if (obj.find("position") != obj.end()) { object->position = vector_to_vec3(obj["position"]); }
	if (obj.find("normal") != obj.end()) { object->normal = vector_to_vec3(obj["normal"]); }

	json material = obj["material"];
	if(material.find("ambient") != material.end()) { object->hasAmbient = true; object->ambient = vector_to_vec3(material["ambient"]); }
	if (material.find("diffuse") != material.end()) {	object->hasDiffuse = true; object->diffuse = vector_to_vec3(material["diffuse"]); }
	if (material.find("specular") != material.end()) { object->hasSpecular = true; object->specular = vector_to_vec3(material["specular"]); }
	if (material.find("transmissive") != material.end()) { object->isTransmissive = true; object->transmissive = vector_to_vec3(material["transmissive"]); }
	if (material.find("reflective") != material.end()) { object->isReflective = true; object->reflective = vector_to_vec3(material["reflective"]); }
	if (material.find("refraction") != material.end()) { object->isRefractive = true; object->refraction = (float)material["refraction"]; }
	if (material.find("shininess") != material.end()) {  object->shininess = (float)material["shininess"]; }
	if (material.find("gloss") != material.end()) { object->gloss = (float)material["gloss"]; object->isGlossy = true; }
	if (material.find("albedo") != material.end()) { object->albedo = (float)material["albedo"]; object->hasAlbedo = true; }
	if (material.find("roughness") != material.end()) { object->roughness = (float)material["roughness"]; object->hasRoughness = true; }

	
	if (obj["type"] == "mesh")
	{
		json &triangles = obj["triangles"];
		int i = 0;

		for (json::iterator triit = triangles.begin(); triit != triangles.end(); ++triit)
		{

			json &tri = *triit;
			point3 a = vector_to_vec3(tri[0]);
			point3 b = vector_to_vec3(tri[1]);
			point3 c = vector_to_vec3(tri[2]);

			object->triangles[i][0] = vector_to_vec3(tri[0]);
			object->triangles[i][1] = vector_to_vec3(tri[1]);
			object->triangles[i][2] = vector_to_vec3(tri[2]);


			//std::cout << "\nTriangle p1: (" << object->triangles[i][0][0] << ","
				//<< object->triangles[i][0][1] << "," << object->triangles[i][0][2] << ")";

			i++;
		}


		object->numTriangles = i;
		calcCenterMass(object);
	}

	
	
	return object;
}

struct Light* createLight(json li)
{
	Light* light = new Light;

	//set type bools
	if (li["type"] == "point") { light->isPoint = true; light->isDirectional = false; light->isAmbient = false; light->isSpot = false; light->isArea = false; }
	else if(li["type"] == "ambient") { light->isAmbient = true; light->isDirectional = false; light->isPoint = false; light->isSpot = false;light->isArea = false;}
	else if (li["type"] == "directional") { light->isDirectional = true; light->isPoint = false; light->isAmbient = false; light->isSpot = false; light->isArea = false;}
	else if (li["type"] == "spot") { light->isSpot = true; light->isDirectional = false; light->isAmbient = false; light->isPoint = false; light->isArea = false;}
	else if (li["type"] == "area") { light->isSpot = false; light->isDirectional = false; light->isAmbient = false; light->isPoint = false; light->isArea = true; }

	light->color = vector_to_vec3(li["color"]);
	if (li.find("position") != li.end()) { light->hasPosition = true; light->position = vector_to_vec3(li["position"]); }
	if (li.find("direction") != li.end()) { light->direction = vector_to_vec3(li["direction"]); }
	if (li.find("cutoff") != li.end()) { light->cutoff = (float)(li["cutoff"]); }
	if (li.find("radius") != li.end()) { light->radius = (float)(li["radius"]); }
	
	return light;
}

void convertScene()
{
	int i = 0;
	for (json::iterator it = sceneObjects.begin(); it != sceneObjects.end(); ++it)
	{
		json object = *it;
		objects[i] = createObject(object);
		if (objects[i]->isMesh) { transformMesh(i); }
		if (objects[i]->isMesh) { objects[i]->bvh = createBoundingVolume(i, 0, objects[i]->numTriangles-1); }
		i++;
	}
	numObjects = i;
	i = 0;
	
	for (json::iterator it = sceneLights.begin(); it != sceneLights.end(); ++it)
	{
		json light = *it;
		lights[i] = createLight(light);
		i++;
	}
	numLights = i;

	for (int j = 0; j < i;j++)
	{
		std::cout << "\nLight color: (" << lights[j]->color[0] << "," << lights[j]->color[1] << "," << lights[j]->color[2] << ")";
	}
}


/****************************************************************************/

void choose_scene(char const *fn) {
	if (fn == NULL) {
		std::cout << "Using default input file " << PATH << "c.json\n";
		fn = "c";
	}

	std::cout << "Loading scene " << fn << std::endl;
	
	std::string fname = PATH + std::string(fn) + ".json";
	std::fstream in(fname);

	if (!in.is_open()) {
		std::cout << "Unable to open scene file " << fname << std::endl;
		exit(EXIT_FAILURE);
	}
	
	in >> scene;


	
	json camera = scene["camera"];
	sceneLights = scene["lights"];
	sceneObjects = scene["objects"];

	convertScene();

	// these are optional parameters (otherwise they default to the values initialized earlier)
	if (camera.find("field") != camera.end()) {
		fov = camera["field"];
		std::cout << "Setting fov to " << fov << " degrees.\n";
	}
	if (camera.find("background") != camera.end()) {
		background_colour = vector_to_vec3(camera["background"]);
		std::cout << "Setting background colour to " << glm::to_string(background_colour) << std::endl;
	}
}

bool inShadow(point3 p, point3 L, point3 &percent)
{
	point3 d = L - p;

	bool completeShadow = false;

	for (int i = 0; i < numObjects && !completeShadow; i++)
	{
		if (objects[i]->isPlane)
		{
			point3 pos = objects[i]->position;
			point3 N = objects[i]->normal;

			float denom = glm::dot(N, L);
			float num = glm::dot(N, (pos - p));
			float t = num / denom;

			if (denom >= tolerance) 
			{	if (t > tolerance) 
				{ 
				if (objects[i]->isTransmissive)
				{
					percent = percent * objects[i]->transmissive;
				}
				else { completeShadow = true; }
					
				} 
			}
		}

		if (objects[i]->isSphere)
		{

			point3 c = objects[i]->position;
			float r = objects[i]->radius;

			float discriminant = (pow(glm::dot(d, (p - c)), 2) - glm::dot(d, d)*(glm::dot((p - c), (p - c)) - pow(r, 2)));

			if (discriminant >= 0)
			{
				float t = (glm::dot(-d, p - c) - sqrt(discriminant)) / glm::dot(d, d);
				if (t >= tolerance) 
				{ 
					if (objects[i]->isTransmissive)
					{
						percent = percent * objects[i]->transmissive;
					}
					else { completeShadow = true; }
					//return true; 
				}
			}
		}

		if (objects[i]->isMesh)
		{
			float tmin = 0;
			int tri = -1;
			traverseBoundingVolumes(p,d,objects[i]->bvh,tmin,tri,false);

			if (tri > -1)//if it hit something
			{
				if (objects[i]->isTransmissive) { percent = percent * objects[i]->transmissive; }
				else { completeShadow = true; }
			}
		}
	}

	return completeShadow;
}



point3 calcReflection(point3 point, point3 e, point3 n, point3 r, int level, bool pick)
{
	point3 finalColour = point3(0, 0, 0);
	if (level < maxReflect)
	{
		point3 d = r;
		int closest = -1;
		point3 N = point3(0, 0, 0);
		float nearestT = nearT(point, d, closest, N, false);
		float t = nearestT;

		if (closest > -1) {
			if (objects[closest]->isPlane)
			{

				if (pick) { std::cout << "\nreflecting plane"; }
				point3 p = (point + t*d);
				finalColour = calcLighting(d, point, p, N, closest, level + 1, pick);

			}

			if (objects[closest]->isSphere)
			{
				point3 c = objects[closest]->position;

				point3 p = point3(point + (d*t));

				point3 N = glm::normalize(p - c);
				if (pick) { std::cout << "\nreflecting sphere"; }
				finalColour = calcLighting(d, point, p, N, closest, level + 1, pick);

			}

			if (objects[closest]->isMesh)
			{

				point3 x = (point + t*d);
				finalColour = calcLighting(d, point, x, N, closest, level + 1, pick);
			}
		}
	}
	else { finalColour = background_colour; }
	return finalColour;
}

point3 calcGlossyReflection(point3 point, point3 e, point3 n, float a, int level, bool pick)
{

	point3 view = glm::normalize(e - point);
	point3 r = glm::normalize(2 * glm::dot(n, view)*n - view);
	point3 k = point3(0,0,0);
	//orthonormal basis
	if (r.x > r.y && r.x > r.z) { k = point3(1, 0, 0); }
	else if (r.y > r.x && r.y > r.z) { k = point3(0, 1, 0); }
	else { k = point3(0, 0, 1); }
	
	point3 u = glm::normalize((glm::cross(k,r)));
	point3 v = glm::normalize((glm::cross(r, u)));


	point3 avg = point3(0, 0, 0);
	//point3 avg = calcReflection(point, e, n, r, level + 1, pick);// calc at least 1 perfect reflection

	for (int i = 0; i < glossSamples;i++)
	{
		float rand0 = (rand() / (float)RAND_MAX);
		float rand1 = (rand() / (float)RAND_MAX);

		//if (pick) { std::cout << "\nrand0 : " << rand0 << "\tRand1 : " << rand1; }

		float uu = -a / 2 + (rand0*a);
		float vv = -a / 2 + (rand1*a);

		point3 rr = glm::normalize((r + (point3(uu*u[0], uu*u[1], uu*u[2]) + point3(vv*v[0], vv*v[1], vv*v[2]))));

		avg += calcReflection(point, e, n, rr, level+1, pick);

	}
	avg = avg / point3((float)glossSamples, (float)glossSamples, (float)glossSamples);
	return avg;
}


point3 calcRefraction(point3 d, point3 e, point3 point, point3 n,int obj, float ni, float nr, int level, bool pick)
{
	


		float vin = glm::dot(d, n);
		point3 vr = point3(0, 0, 0);
		point3 finalColour = point3(0, 0, 0);

		if (pow(vin, 2) < 1 - pow(ni, 2) / pow(nr, 2))
		{
			vr = ni*(d - n*(vin)) / nr - n * (sqrt(1 - pow(ni, 2)*(1 - pow(vin, 2)) / pow(nr, 2)));
			vr = glm::normalize(vr);
			d = vr;

			if (pick)
			{
				std::cout << "\nVR: (" << vr[0] << "," << vr[1] << "," << vr[2] << ")";
				//std::cout << "\nt: " << t;
			}
		

		int closest = -1;
		point3 N = point3(0, 0, 0);
		float nearestT = nearT(point, d, closest, N, false);
		float t = nearestT;

		

		if (closest > -1)
		{
			
				if (objects[closest]->isPlane)
				{

					if (pick) { std::cout << "\nrefracting plane"; }
					point3 p = (point + t*d);
					finalColour = calcLighting(d, point, p, N, closest, level, pick);

				}

				if (objects[closest]->isSphere)
				{

					point3 c = objects[closest]->position;

					point3 p = point3(point + (d*t));

					if (closest == obj) { finalColour = calcRefraction(d, point, p, -N, obj, level, nr, ni, pick); }
					else {

						point3 N = glm::normalize(p - c);
						if (pick) { std::cout << "\nrefracting sphere"; }
						finalColour = calcLighting(d, point, p, N, closest, level, pick);
					}

				}


				if (objects[closest]->isMesh)
				{
					
						point3 x = (point + t*d);
						if (closest == obj) { finalColour = calcRefraction(d, point, x, -N, obj, level, nr, ni, pick); }
						else {
						finalColour = calcLighting(d, point, x, N, closest, level, pick);
					}
				}
		}
		else { finalColour = background_colour; }
		}
		else { return background_colour; }

return finalColour;
}


point3 calcTransmission(point3 d, point3 e, point3 point, point3 N, int obj, int level, bool pick)
{
	int closest = -1;
	float t = nearT(point, d, closest, N, pick);
	
	point3 finalColour = point3(0, 0, 0);
	if (closest > -1)
	{
		if (objects[closest]->isPlane)
		{
			point3 p = (point + t*d);
			finalColour = calcLighting(d, point, p, N, closest, level, pick);
		}

		if (objects[closest]->isSphere) {

			point3 c = objects[closest]->position;

			point3 p = point3(point + (d*t));

			point3 N = glm::normalize(p - c);

			finalColour = calcLighting(d, point, p, N, closest, level, pick);
		}

		if (objects[closest]->isMesh)
		{
			point3 x = (point + t*d);
			finalColour = calcLighting(d, point, x, N, closest, level, pick);

		}
	}
	return finalColour;
}
point3 calcAreaLightPercentage(point3 p, int light, bool pick)
{
	float rad = lights[light]->radius;
	point3 avg = point3(0, 0, 0);
	point3 finalPercentage = point3(1, 1, 1);
	point3 pos = lights[light]->position;

	for (int i = 0; i < areaSamples;i++)
	{
		float x = (rand() / (float)RAND_MAX * 2 * rad) - rad;
		float y = (rand() / (float)RAND_MAX * 2 * rad) - rad;
		float z = (rand() / (float)RAND_MAX * 2 * rad) - rad;

		point3 randPoint = point3(x, y, z);
		randPoint = pos + randPoint;
		randPoint = glm::normalize(randPoint);
		if (pick) { std::cout << "\nRand : (" << randPoint[0] << "," << randPoint[1] << "," << randPoint[2] << ")"; }
		point3 L = glm::normalize(randPoint - p);

		point3 per = point3(1, 1, 1);

		if (!inShadow(p, randPoint, per))
		{
			avg += point3(1, 1, 1);//*per; // * 1 if there is no transmissive object
		}
	}
	avg = avg / point3((float)areaSamples, (float)areaSamples, (float)areaSamples);
	if (pick) { std::cout << "\nAvg : (" << avg[0] << "," << avg[1] << "," << avg[2] << ")"; }
	return avg;
}


point3 calcLight(point3 d, point3 e, point3 p, point3 N, int obj, int light, int level, bool pick) 
{

	point3 finalColour = point3(0, 0, 0);
	point3 percent = point3(1, 1, 1);
	point3 col = lights[light]->color;
	point3 L = point3(0, 0, 0);
	point3 pos = point3(0, 0, 0);
	point3 dir = point3(0, 0, 0);

	pos = lights[light]->position;

	float cutoff = 0;
	float spotL = 0;

	if (lights[light]->isDirectional) { dir = lights[light]->direction; L = glm::normalize(-dir); }
	if (lights[light]->isPoint || lights[light]->isArea) { L = glm::normalize(pos - p); }
	if (lights[light]->isSpot)
	{
		dir = lights[light]->direction;
		L = glm::normalize(pos - p);
		if (pick)
		{
			std::cout << "\npos : (" << pos[0] << "," << pos[1] << "," << pos[2] << ")";
			std::cout << "\np : (" << p[0] << "," << p[1] << "," << p[2] << ")";
			std::cout << "\nL : (" << L[0] << "," << L[1] << "," << L[2] << ")";
		}
		cutoff = lights[light]->cutoff;
		cutoff = glm::radians(cutoff);
		spotL = glm::dot(dir, -L);
		spotL = acos((spotL));

		if (pick)
		{
			std::cout << "\nspotL :" << spotL;
			std::cout << "\ncutoff :" << cutoff;
		}
	}

	if (!inShadow(p, L, percent) || lights[light]->isArea)
	{
		
		if (lights[light]->isDirectional || lights[light]->isSpot || lights[light]->isPoint || lights[light]->isArea)
		{

			float NL = (glm::dot(N, L));
			if (NL < 0) { NL = 0; }

			if (!lights[light]->isSpot || (lights[light]->isSpot && spotL < cutoff))
			{
				if (objects[obj]->hasDiffuse)
				{
					point3 diffuse = objects[obj]->diffuse;
					point3 diff = (diffuse*col)*(NL);
				
					if (objects[obj]->hasAlbedo && objects[obj]->hasRoughness)
					{
						//oren-nayar
						float roughness = pow(objects[obj]->roughness, 2);
						float albedo = objects[obj]->albedo;

						point3 view = glm::normalize(e - p);

						float NV = glm::dot(N, view);

						float alpha = acos(std::max(NL, NV));
						float beta = acos(std::min(NL, NV));

						float gamma = glm::dot(view - N * NV, L - N * NL);

						float C1 = 1.0f - 0.5f*(roughness / roughness + 0.33f);
						float C2 = 0.45f * (roughness / roughness + 0.09f);

						if (gamma >= 0) { C2 *= sin(alpha); }
						else { C2 *= (sin(alpha) - pow((2.0f*beta) / glm::pi<float>(), 3)); }

						float power = (4.0f * alpha * beta) / (glm::pi<float>() * glm::pi<float>());
						float C3 = 0.125f * (roughness + 0.09f) * power * power;

						float A = gamma * C2 * tan(beta);
						float B = (1.0f - fabs(gamma)) * C3 * tan((alpha + beta) / 2);

						float L1 = std::max(0.0f, NL) * (C1 + A + B);

						float betaPI = 2.0f * beta / glm::pi<float>();

						float L2 = 0.17f * NL * (roughness / (roughness + 0.13f)) * (1.0f - gamma * betaPI * betaPI);

						point3 oren = diffuse * (L1 + L2);
						finalColour += oren*NL*col;

					}
					else
					{
						
						finalColour += diff;
					}
				}

				if (objects[obj]->hasSpecular)
				{

					if (NL > 0.0)
					{
						point3 r = glm::normalize(2 * (NL)*N - L);

						point3 v = glm::normalize(e - p);

						float RV = glm::dot(r, v);
						if (RV < 0) { RV = 0; }

						point3 specular = objects[obj]->specular;
						float shine = objects[obj]->shininess;

						if (objects[obj]->hasRoughness)
						{
							//ward anisotropic
							float aX = objects[obj]->roughness;
							float aY = aX;

							point3 V = glm::normalize(e - p);
							point3 H = glm::normalize((L) + (V));
							float NV = std::max((float)0.0, glm::dot(N, V));

							if (NV > 0) {

								point3 T = point3(0, 0, 0);
								point3 B = point3(0, 0, 0);
								findOrthonormalVectors(N, T, B);

								float HT = std::max((float)0.0, glm::dot(H, T));
								float HB = std::max((float)0.0, glm::dot(H, B));
								float HN = std::max((float)0.0, glm::dot(H, N));

								float beta = -2.0 *(pow(HT / aX, 2.0) + pow(HB / aY, 2.0)) / (1.0 + HN);
								float denom = std::max(FLT_MIN, (float)(4.0 * glm::pi<float>() * aX * aY * sqrt(NL * NV)));

								point3 spec = objects[obj]->specular;

								spec[0] = (spec[0] * NL / denom) * exp(beta);
								spec[1] = (spec[1] * NL / denom) * exp(beta);
								spec[2] = (spec[2] * NL / denom) * exp(beta);

								point3 ward = spec;

								ward = glm::clamp(ward, point3(0, 0, 0), point3(1, 1, 1));

								if (pick) { std::cout << "\nWard: (" << ward[0] << "," << ward[1] << "," << ward[2] << ")"; }

								specular = ward;
							}
						}



						if (objects[obj]->hasRoughness) {
							finalColour += specular*col;//*pow(RV, shine);
						}
						else { finalColour += (specular*col)*pow(RV, shine); }

					}
				}
			}

			finalColour = finalColour*percent;
		}

		if (lights[light]->isArea) 
		{ 
			finalColour = finalColour *calcAreaLightPercentage(p, light, pick);
		}
	} 

	if (lights[light]->isAmbient && objects[obj]->hasAmbient)
	{
		point3 ambient = objects[obj]->ambient;
		finalColour += (ambient*col);
	}
	
	return finalColour;
}

point3 calcLighting( point3 d, point3 e, point3 p, point3 N, int obj,int level, bool pick)
{

	point3 finalColour = point3(0, 0, 0);

	for(int i = 0; i < numLights;i++)
	{
		finalColour += calcLight(d, e, p, N, obj, i, level, pick);
	}

	if(objects[obj]->isReflective)
	{
		point3 reflective = objects[obj]->reflective;
		if (objects[obj]->isGlossy) 
		{
			float gloss = objects[obj]->gloss;
			if (level < maxReflect) { finalColour += calcGlossyReflection(p, e, N, gloss, level, pick)*reflective; }
		}
		else {
			
			point3 v = e - p;
			point3 r = 2 * glm::dot(N, v)*N - v;
			if (level < maxReflect) { finalColour += (calcReflection(p, e, N, r, level, pick)*reflective); }
			if (pick)
			{
				std::cout << "\nreflective: (" << reflective[0] << "," << reflective[1] << "," << reflective[2] << ")";
				std::cout << "\nfc: (" << finalColour[0] << "," << finalColour[1] << "," << finalColour[2] << ")";
			}
		}
	}

	if(objects[obj]->isTransmissive)
	{
		point3 trans = objects[obj]->transmissive;
		if(objects[obj]->isRefractive)
		{
			float refraction = objects[obj]->refraction;
			point3 refractedColour = calcRefraction(d, e, p, N,obj,level, refraction,AIR, pick);//AIR THEN REFRACTION
			finalColour = (point3(1, 1, 1) - trans)*finalColour + trans*refractedColour;
		}
		else {

			point3 transColour = calcTransmission(d, e, p, N, obj, level, pick);
			finalColour = (point3(1, 1, 1) - trans)*finalColour + trans*transColour;

			if (pick)
			{

				//std::cout << "\nTransmissive: (" << trans[0] << "," << trans[1] << "," << trans[2] << ")";
				//std::cout << "\nfc: (" << finalColour[0] << "," << finalColour[1] << "," << finalColour[2] << ")";
			}

		}
	}

	return finalColour;
}

bool traverseBoundingVolumes(point3 e, point3 d, BVHnode* bvh, float &tmin, int& tri, bool pick)
{


	int mesh = bvh->mesh;
	
	{
		//kajiya
		
		float tnear = -100;
		float tfar = 100;

		for (int i = 0; i < 3; i++)
		{
			if (d[i] == 0) 
			{
				if (e[i] < bvh->bvMin[i] || e[i] > bvh->bvMax[i])
				{
					return false;
				}
			}
			else
			{
				float t1 = (bvh->bvMin[i] - e[i]) / d[i];
				float t2 = (bvh->bvMax[i] - e[i]) / d[i];

				if (t1 > t2) { std::swap(t1, t2); }
				if (t1 > tnear) { tnear = t1; }
				if (t2 < tfar) { tfar = t2; }


				if (tnear > tfar) { return false; }
				if (tfar < 0) { return false; }
			}
		}

		tmin = tnear;

		if (bvh->isLeaf)
		{
			tri = bvh->rangeMax;
			if (pick) {
				std::cout << "\nhit leaf! Tri: " << tri; 
			}

		}
		else {

			float tminleft = 100;
			int trileft = -1;
			float tleft = 0;
			bool lefthastri = false;

			traverseBoundingVolumes(e, d, bvh->left, tminleft, trileft, pick);

			if (trileft > -1)
			{

				point3 a = objects[mesh]->triangles[trileft][0];
				point3 b = objects[mesh]->triangles[trileft][1];
				point3 c = objects[mesh]->triangles[trileft][2];

				point3 N = glm::normalize(glm::cross(b - a, c - a));

				tleft = glm::dot(N, (a - e)) / glm::dot(N, d);

				point3 x = (e + tleft*d);

				float x0 = glm::dot(glm::cross((b - a), (x - a)), N);
				float x1 = glm::dot(glm::cross((c - b), (x - b)), N);
				float x2 = glm::dot(glm::cross((a - c), (x - c)), N);

				if (x0 > 0 && x1 > 0 && x2 > 0) 
				{ 
					if (tleft > tolerance)
					{
						tri = trileft;
						lefthastri = true;
					}
				}
			}
				
			float tminright = 100;
			int triright = -1;
			float tright = 0;
			bool righthastri = false;

			traverseBoundingVolumes(e, d, bvh->right, tminright, triright, pick);
			if (triright > -1 )
			{
				point3 a = objects[mesh]->triangles[triright][0];
				point3 b = objects[mesh]->triangles[triright][1];
				point3 c = objects[mesh]->triangles[triright][2];

				point3 N = glm::normalize(glm::cross(b - a, c - a));

				tright = glm::dot(N, (a - e)) / glm::dot(N, d);
				point3 x = (e + tright*d);

				float x0 = glm::dot(glm::cross((b - a), (x - a)), N);
				float x1 = glm::dot(glm::cross((c - b), (x - b)), N);
				float x2 = glm::dot(glm::cross((a - c), (x - c)), N);

				if (x0 > 0 && x1 > 0 && x2 > 0)
				{
					if (tright > tolerance)
					{
						tri = triright;
						righthastri = true;
					}
				}
			}
			
			if (lefthastri && righthastri) // if they both hit a tri, which one is closer
			{
				if (tleft < tright) { tri = trileft; }
				else if (tright < tleft) { tri = triright; }
			}
		}
		return true;
	}
}

//finds nearest t value
float nearT(point3 e, point3 d, int &obj, point3 &normal, bool pick)
{
	float nearT = FLT_MAX;

	for(int i = 0; i < numObjects;i++)
	{
		if(objects[i]->isPlane)
		{

			point3 pos = objects[i]->position;
			point3 N = objects[i]->normal;

			float t = glm::dot(N, (pos - e)) / glm::dot(N, d);
			if (pick) { std::cout << "\nplaneT: " << t; }
			if (t > tolerance) { if (t < nearT) { nearT = t; obj = i; normal = N; } }
			//if (pick) { std::cout << "\nType: " << obj["type"]; }
		}

		if(objects[i]->isSphere)
		{
			point3 c = objects[i]->position;
			float r = objects[i]->radius;

			float discriminant = (pow(glm::dot(d, (e - c)), 2) - glm::dot(d, d)*(glm::dot((e - c), (e - c)) - pow(r, 2)));
			float t2 = (glm::dot(-d, e - c) + sqrt(discriminant)) / glm::dot(d, d);
			float t = (glm::dot(-d, e - c) - sqrt(discriminant)) / glm::dot(d, d);//??
			if (pick) { std::cout << "\nspheret: " << t; }
			if (pick) { std::cout << "\nspheret2: " << t2; }
			if (discriminant >= 0) {
				if (t >= tolerance) { if (t < nearT) { nearT = t; obj = i;} }
				else { if (t2 >= tolerance) { if (t2 < nearT) { nearT = t2; obj = i;} } }
			}
		}

		if(objects[i]->isMesh)
		{
			int tri = -1;
			float tmin = 0;
			traverseBoundingVolumes(e, d, objects[i]->bvh, tmin, tri, pick);

			if (pick)
			{
				tmin = tmin;
			}
			if (tri > -1)
			{
			
				point3 a = objects[i]->triangles[tri][0];
				point3 b = objects[i]->triangles[tri][1];
				point3 c = objects[i]->triangles[tri][2];

				//std::cout << "\nTriangle a: (" << object[i]->triangles[i][0][0] << ","
				//<< object->triangles[i][0][1] << "," << object->triangles[i][0][2] << ")";
				point3 N = glm::normalize(glm::cross(b - a, c - a));

				float t = glm::dot(N, (a - e)) / glm::dot(N, d);
				if (pick) { std::cout << "\nmeshT: " << t; }

				point3 x = (e + t*d);

				float x0 = glm::dot(glm::cross((b - a), (x - a)), N);
				float x1 = glm::dot(glm::cross((c - b), (x - b)), N);
				float x2 = glm::dot(glm::cross((a - c), (x - c)), N);

				if (x0 > 0 && x1 > 0 && x2 > 0)
				{
					if (t > tolerance) { if (t < nearT) { nearT = t; obj = i; normal = N; } } //N ****
				}
				
			}
		}
	}
	if (pick) { std::cout << "\nnearestT: " << nearT; }
	return nearT;
}

bool trace(const point3 &e, const point3 &s, colour3 &colour, bool pick) {
	// NOTE 1: This is a demo, not ray tracing code! You will need to replace all of this with your own code...
  // NOTE 2: You can work with JSON sceneObjects directly (like this sample code), read the JSON sceneObjects into your own data structures once and render from those (probably in choose_scene), or hard-code the sceneObjects in your own data structures and choose them by name in choose_scene; e.g. choose_scene('e') would pick the same scene as the one in "e.json". Your choice.
  // If you want to use this JSON library, https://github.com/nlohmann/json for more information. The code below gives examples of everything you should need: getting named values, iterating over arrays, and converting types.

	point3 d = (s - e);
	int closest = -1;
	point3 N = point3(0, 0, 0);
	float nearestT = nearT(e,d,closest,N,pick);
	
	float t = nearestT;
	int level = 0;

	if (closest > -1)
	{
		if (objects[closest]->isPlane)
		{

			point3 pos = objects[closest]->position;
			point3 N = objects[closest]->normal;

			point3 p = (e + t*d);
			colour = calcLighting(d, e, p, N, closest, level, pick);
			return true;

		}

		if (objects[closest]->isSphere) {

			point3 c = objects[closest]->position;

			point3 p = point3(e + (d*t));

			point3 N = glm::normalize(p - c);

			colour = calcLighting(d, e, p, N, closest, level, pick);

			return true;

		}
		if (objects[closest]->isMesh)
		{

			point3 x = (e + t*d);
			///////////////XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX//////////////////
			colour = calcLighting(d, e, x, N, closest, level, pick);
			return true;
		}
	}
	else { colour = background_colour; }
	return false;
}