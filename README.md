# Framework for surface mesh processing （C++Framework for  Assignments of GAMES301 ）

**The program , which runs in the Qt application, is used to load meshes and render using OpenGL.**


## Dependent External Libraries
* [Qt](https://www.qt.io/), Recommended version: 5.13.0
## Alternative External Libraries
* [Eigen](http://eigen.tuxfamily.org/)
* [OpenMesh](https://www.openmesh.org/), Recommended version: the latest 8.1(at Oct. 2020)

## Usage


```
cd Surface_Framework_Cmake
mkdir build && cd build
cmake -A x64 ..
```


Open **SurfaceFrameworkCmake.sln**, select **SurfaceFrameworkCmake** as launch project, and run.

## Supported File Formats

.obj .off .ply .stl



# Note

- 使用vspkg导入Eigen库

## TODO

- 修改固定管线为可编程管线
- 纹理贴图
- 将参数化结果映射回输入模型
- Hw2中最小步长过小
- SD以方形开始时出现nan







