# A MeshLib demo

This C++ project framework is used to help students to debug geometric algorithms. It contains a halfedge data structure library `MeshLib` and an simple opengl viewer.

## System

The code is only tested on Windows, but it should work on Linux and Mac with minor midifications. If there is any problem on the latter two platforms, please let me know.

## Dependencies

1. `MeshLib`a mesh library based on halfedge data structure.
2. `freeglut`, a free-software/open-source alternative to the OpenGL Utility Toolkit (GLUT) library.

## Directory Structure

``` txt
3rdparty       -- The `MeshLib` and `freeglut` library.
demo_project   -- The source code of `demo_project`.
data           -- Some models and texture images.
resources      -- Some pictures needed by README.md.
CMakeLists.txt -- CMake configuration file.
```

## Configuration

### Windows

1. Install `MeshLib`, `freeglut`, they are included in `3rdparty`.

> Note: You need to add `YOUR_PATH/freeglut/bin/x64` to the system or user `PATH` variable in order to ensure other executable programs can find `freeglut.dll`.

2. Install [CMake](https://cmake.org/download/).

3. Download the source code of the C++ framework.
> E.x. I create a folder `projects` in `C:/`, then unzip the code there.

4. Configure and generate the project for Visual Studio.
  
> ``` bash
> cd MyProjects
> mkdir build
> cd build
> cmake ..
> ```
> *One can also finish this step using CMake GUI.*

5. Open the \*.sln using Visual Studio, and complie `DemoProject` and `INSTALL`.

6. Run the executable program.
> E.x. 
> ``` bash
> cd ../bin/
> ./DemoProject.exe ../data/kitten.m ../data/texture_checker.bmp
> ./DemoProject.exe ../data/Alex_rgb.m
> ```
> If you can see the following results, then it means that the configuration succeeded. 
> 
> ![kitten](kitten.png) ![Alex](Alex.png)
