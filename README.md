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

```c++
//cpp file
#include <Utils.h>
namespace str=lyxutils::str_utils;//simplify namespaces
```

