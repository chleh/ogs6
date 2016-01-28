function(documentationProjectFilePutIntoPlace p)
    file(RELATIVE_PATH relative_path ${DocumentationProjectFileInputDir} ${p})
    get_filename_component(dir_name ${relative_path} DIRECTORY)

    get_filename_component(otagname ${relative_path} NAME_WE)
    if (otagname MATCHES ^__)
        # if the file name starts with an underscore, then this files is
        # the "table of contents of the current directory
        
        file(MAKE_DIRECTORY "${DocumentationProjectFileBuildDir}/${dir_name}")

        set(postfix "# Parameters used at this level\n\n")

        # gather other parameter files
        # the loop below will effects a page hierarchy to be built
        file(GLOB param_files ${DocumentationProjectFileInputDir}/${dir_name}/*)
        list(SORT param_files)
        foreach(pf ${param_files})
            get_filename_component(rel_pf ${pf} NAME_WE)

            # if the file name starts with an underscore, then this
            # is the "table of contents" file already processed outside
            # of this loop
            if (NOT rel_pf MATCHES ^__)
                if(IS_DIRECTORY ${pf})
                    set(pf_tagname ${rel_pf})
                else()
                    string(SUBSTRING ${rel_pf} 2 -1 pf_tagname)
                endif()

                if ("${dir_name}" STREQUAL "") # toplevel dir must be trested slightly different
                    set(pf_tagpath "${pf_tagname}")
                else()
                    set(pf_tagpath "${dir_name}/${pf_tagname}")
                    string(REPLACE "/" "__" pf_tagpath "${pf_tagpath}")
                endif()
                message("  t.o.c. entry ${pf_tagpath}")

                set(postfix "${postfix} - \\subpage ogs_project_file_parameter__${pf_tagpath}\n")

                if (NOT IS_DIRECTORY "${pf}")
                    documentationProjectFilePutIntoPlace("${pf}")
                endif()
            endif()
        endforeach()
    else()
        set(postfix "")
    endif()

    string(SUBSTRING ${otagname} 2 -1 tagname)
    if (dir_name STREQUAL "") # toplevel dir must be trested slightly different
        set(tagpath "${tagname}")
    else()
        if (otagname MATCHES ^__) # treat "table of contents" file special
            string(REPLACE "/" "__" tagpath "${dir_name}")
        else()
            string(REPLACE "/" "__" tagpath "${dir_name}/${tagname}")
        endif()
    endif()
    message("  child param  ${tagpath}")

    set(pagetitle "[tag]&emsp;${tagname}")

    # read, augment, write file content
    file(READ ${p} content)
    set(content "/*! \\page ogs_project_file_parameter__${tagpath} ${pagetitle}\n${content}\n\n${postfix}\n*/\n")
    string(REGEX REPLACE .md$ .dox output_file "${DocumentationProjectFileBuildDir}/${relative_path}")
    file(WRITE "${output_file}" "${content}")
endfunction()


set(DocumentationProjectFileBuildDir ${PROJECT_BINARY_DIR}/DocProjectFileAux)
set(DocumentationProjectFileInputDir ${PROJECT_SOURCE_DIR}/Documentation/ProjectFile)

# remove old output
if (IS_DIRECTORY ${DocumentationProjectFileBuildDir})
	file(REMOVE_RECURSE ${DocumentationProjectFileBuildDir})
endif()

file(GLOB_RECURSE input_paths ${DocumentationProjectFileInputDir}/__*)

foreach(p ${input_paths})
	message("directory index file ${p}")
	documentationProjectFilePutIntoPlace(${p})
endforeach()
