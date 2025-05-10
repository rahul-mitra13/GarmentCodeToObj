# Dataset and Purpose: 

We use the dataset [Generating Datasets of 3D Garments with Sewing Patterns](https://zenodo.org/records/5267549) from [Maria Korosteleva](https://korosteleva.com/) and [Sung-Hee Lee](https://lava.kaist.ac.kr/?page_id=41).

This code converts json specifications in this dataset to a trianglutated mesh with appropriate "stitched" together vertices at panel borders. The mesh is stored in an obj file - `mesh.obj`. This allows for downstream applications/processing that may require a mesh setting. 

# Usage 

All the requirements are listed in `requirements.txt`. Here are the build instruction in Linux/MacOS systems. I recommend making a virtual enviroment. 

```
python3 -m venv venv
source venv/bin/activate
pip3 install -r requirements.txt
python3 garmentJsonToObj.py [GarmentCode JSON file] [Max triangle area constraint]
```

The JSON file is named `specification.json` in their dataset and the max triangle area is the constraint imposed by the [Triangle](https://www.cs.cmu.edu/~quake/triangle.switch.html) library (what we use for a Constrained Delaunay Triangulation of the panels).

# Output 

[Polyscope](https://polyscope.run/) is used to render the output. For example, here is the output from using `specification.json` from `jumpsuit_sleeveless_1KH03EYORI` model in their dataset (max triangle area was set to 1.0).

<img width="731" alt="Screenshot 2025-05-09 at 3 12 04 PM" src="https://github.com/user-attachments/assets/dc30cf65-c07e-4106-86c9-e861ec3d9d0e" />



An obj file is also saved, `mesh.obj` which, in addition to the standard `v` and `f` labels, also stores an `m` label which stores pairs of "stitched together" vertices. 
