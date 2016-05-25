if(DOXYGEN_FOUND)

    option(DOCS_GENERATE_DOCSET "Generate Dash Docsets." OFF)

    set(DOT_FOUND "NO")
    if(DOXYGEN_DOT_FOUND)
        set(DOT_FOUND "YES")
    endif()

    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Generating source code documentation with Doxygen." VERBATIM)

    # Defaults
    set(DOCS_GENERATE_TREEVIEW_STRING "YES" CACHE INTERNAL "")
    set(DOCS_DISABLE_INDEX_STRING "NO" CACHE INTERNAL "")
    set(DOCS_GENERATE_DOCSET_STRING "NO" CACHE INTERNAL "")
    set(DOCS_SEARCHENGINE_STRING "YES" CACHE INTERNAL "")

    # Dash Docsets
    if(DOCS_GENERATE_DOCSET)
        find_program(DOCSETUTIL_TOOLPATH docsetutil
            PATH /Applications/Xcode.app/Contents/Developer/usr/bin)
        if(NOT DOCSETUTIL_TOOLPATH)
            message(FATAL_ERROR "docsetutil required for Docset-generation!")
        endif()
        set(DOCS_GENERATE_TREEVIEW_STRING "NO" CACHE INTERNAL "")
        set(DOCS_DISABLE_INDEX_STRING "YES" CACHE INTERNAL "")
        set(DOCS_GENERATE_DOCSET_STRING "YES" CACHE INTERNAL "")
        set(DOCS_SEARCHENGINE_STRING "NO" CACHE INTERNAL "")
        add_custom_command(TARGET doc POST_BUILD
            COMMAND make
            COMMAND mv org.doxygen.Project.docset ogs6.docset
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/docs
            COMMENT "Generating docset ...")
        configure_file(Documentation/DocsetFeed.xml.in ${PROJECT_BINARY_DIR}/docs/ogs6.xml)
    endif()

    configure_file(Documentation/Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile)

	add_custom_target(internal_pre_doc
		${CMAKE_COMMAND}
		-DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
		-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
		-P ${PROJECT_SOURCE_DIR}/scripts/cmake/DocumentationProjectFile.cmake
		WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
		COMMENT "Generating project file documentation hierarchy." VERBATIM)
	add_dependencies(doc internal_pre_doc)
endif()
