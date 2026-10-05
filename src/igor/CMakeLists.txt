### CMakeLists.txt ---

# Layers in dependency order (doc/LAYER_REFACTORING_PROPOSAL.md, §2):
# Core and Math depend on nothing; Alignment on Core; Model on Core and Math;
# Inference and Generation on Model; Streaming on Core.
add_subdirectory(Core)
add_subdirectory(Math)
add_subdirectory(Alignment)
add_subdirectory(Model)
add_subdirectory(Generation)
add_subdirectory(Inference)
add_subdirectory(Streaming)

######################################################################
### CMakeLists.txt ends here
