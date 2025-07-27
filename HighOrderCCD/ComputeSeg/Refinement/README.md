# README

## 包含进VS工程

添加AVX_utils.h和CPS_xxx.h文件到VS工程中。然后设置预定义项

* 使用普通CPU指令集，需要依赖CGAL，在VS工程设置里添加预定义项 USE_VEC
* **（推荐）**使用AVX指令集，不依赖CGAL，在VS工程设置里添加预定义项 USE_AVX

添加预定义成功与否，可以看 CPS_Vector.h中哪一段代码被激活。

在需要使用的地方，包含"CPS_AABBTree.h"即可，注意有命名空间ClosestPointSearch。

## 接口

提供了接受OpenMesh的Mesh的构造函数。

最近点查询函数接受OpenMesh的Vec3d，可返回最近点、最近点距离、最近的三角形。详见AABBTree类的public下，
以closest开头的函数。
