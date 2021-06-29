LibraryManager_Add(${PROJECT_NAME} NAMESPACE molpro)
add_subdirectory(molpro)
LibraryManager_Install(${PROJECT_NAME} PKG_CONFIG)
LibraryManager_Export(${PROJECT_NAME} NAMESPACE molpro)

add_executable(main main.cpp)
set_target_properties(main PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
target_link_libraries(main PRIVATE molpro::${PROJECT_NAME})

