add_executable(testStr testStr.cpp)
target_link_libraries(testStr utils)

add_executable(cltest CLTest.cpp)
target_link_libraries(cltest utils)