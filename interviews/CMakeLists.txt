project(interviews)
cmake_minimum_required(VERSION 2.8)

add_executable(t1 test1.cpp test1.h)
add_executable(stl stl.cpp stl.h)
add_executable(ds dsmain.cpp bst.cpp bst.h heap.cpp heap.h trie.cpp trie.h)
add_executable(gl glmain.h glmain.cpp height_union.h height_union.cpp)