# Utils

自己用的一些工具函数。编译安装：

```cmd
mkdir build && cd build
cmake ..
make -j4
sudo make install
```

使用：

```cmake
#CMakeLists.txt
include_directories(/usr/local/include/lyxutils/)
link_directories(/usr/local/lib/)
target_link_libraries(<executable> utils)
```

或者按照一般的cmake查找库的方法使用（目前还不是很会配置cmake的相关config文件，依赖的Eigen库还得在使用的时候再次查找，有待改进）：

```cmake
#CMakeLists.txt
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})
find_package(Utils REQUIRED)
include_directories(${Utils_INCLUDE_DIRS})
link_libraries(${Utils_LIBS})
```

源文件中为方便使用可重定义名字空间为更简短的名字：

```c++
//cpp file
#include <Utils.h>
namespace str=lyxutils::str_utils;//simplify namespaces
```

