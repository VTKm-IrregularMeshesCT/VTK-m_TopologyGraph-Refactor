### Dependencies
vtkm 1.5
vtk 9.0


### Example compilation
cmake -DCMAKE_INSTALL_PREFIX=/home/sc17dd/local-installs/vtk-m-git-refactoring ../
make -j20 install

### Example run
cd build
./cv1k -f ../testing/hydrogen_atom.vtk -o output -t 4
./ContourVisualiser -f ../gridded-04/hh96.txt -o output -t 11 --decompositionType volume
./ContourVisualiser -f ../gridded-04/hh24-space-E.txt -o output -t 15 --decompositionType pactbd

### Evaluation
You can test using the testing/hydrogen_atom.vtk file.
The output files will be saved in build/outputcontour.i.vtk. Compare those against testing/outputcontour.i.vtk.
What you expect is 4 surfaces - two big spheres, a donut shape betwen them, with a small sphere in the middle the donut hole.
See testing/screenshot.png for how it looks like in paraview.




