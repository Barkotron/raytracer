
## COMP 4490 Ray tracer project
#### Jackson Barker


Written and tested using Visual Studio on Windows 10.
To run, type name of scene in project->properties->debugging->command arguments
without '.json'

## Run Times:

gloss (10 samples, 0.9 gloss): ~5:25
gloss (10 samples, 0.2 gloss): ~30 seconds
gloss (24 samples, 0.2 gloss): ~5:14
gloss (16 samples, 0.5 gloss): ~7 minutes
gloss (16 samples, 0.1 gloss): ~40 seconds
area (64 samples): ~5:45
ward: 16 seconds
oren: 16 seconds
i (without BVH): 3:43
i (with BVH): 17 seconds
transform (all): 16 seconds

all gloss screenshots use transform0.
gloss screenshots are named glossXsamplesYglossamount
ex: gloss10samples9 = 24 samples with gloss of 0.9

scene area: scene c with lights replaced by 'area' type lights.
scene transform: scene c with transformations applied (spot light replaced with point light.)
scene oren: left has oren, right has phong - neither has a specular component.
scene ward: left has ward, right has phong.

screenshot naming for transfomations follow the naming below.

### Transformations

transform0:
rotation angle: 90 degrees
rotation vector: (0,1,0.5)
translation: (1,0,0)
scale: (0.5,0.7,0.5)

transform1:
rotation angle: 45 degrees
rotation vector: (0,1,0.5)
translation: (0.5,-0.1,0)
scale: (0.4,0.4,0.4)

transform2:
rotation angle: 110 degrees
rotation vector: (0,0,1)
translation: (0,-0.2,0)
scale: (1,1,1)

transform3:
rotation angle: 25 degrees
rotation vector: (0,0,1)
translation: (0.2,0.0,-10.0)
scale: (5,5,5)
