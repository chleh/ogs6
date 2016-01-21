function(documentationProjectFilePutIntoPlace p)
    file(RELATIVE_PATH relative_path ${DocumentationProjectFileInputDir} ${p})
    get_filename_component(dir_name ${relative_path} DIRECTORY)

    get_filename_component(tagname ${relative_path} NAME_WE)
    if (tagname MATCHES ^_)
        # if the file name starts with an underscore, then this files is
        # the "table of contents of the current directory
        
        file(MAKE_DIRECTORY "${DocumentationProjectFileBuildDir}/${dir_name}")

        set(otagname "${tagname}")
        string(SUBSTRING ${tagname} 1 -1 tagname) # remove _ from start

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
            if (NOT rel_pf MATCHES ^_)
                set(pf_tagname ${rel_pf})

                if ("${dir_name}" STREQUAL "") # toplevel dir must be trested slightly different
                    set(pf_tagpath "${pf_tagname}")
                else()
                    set(pf_tagpath "${dir_name}/${pf_tagname}")
                    string(REPLACE "/" "__" pf_tagpath "${pf_tagpath}")
                endif()

                set(postfix "${postfix} - \\subpage ogs_project_file_parameter__${pf_tagpath}\n")

                if (NOT IS_DIRECTORY "${pf}")
                    documentationProjectFilePutIntoPlace("${pf}")
                endif()
            endif()
        endforeach()
    else()
        set(postfix "")
        set(otagname "")
    endif()

    if (dir_name STREQUAL "") # toplevel dir must be trested slightly different
        set(tagpath "${tagname}")
    else()
        if (otagname MATCHES ^_) # treat "table of contents" file special
            string(REPLACE "/" "__" tagpath "${dir_name}")
        else()
            string(REPLACE "/" "__" tagpath "${dir_name}/${tagname}")
        endif()
    endif()
    message("tp ${tagpath}")

    # read, augment, write file content
    file(READ ${p} content)
    set(content "/*! \\page ogs_project_file_parameter__${tagpath} ${tagname}\n${content}\n\n${postfix}\n*/\n")
    string(REGEX REPLACE .md$ .dox output_file "${DocumentationProjectFileBuildDir}/${relative_path}")
    file(WRITE "${output_file}" "${content}")
endfunction()


set(DocumentationProjectFileBuildDir ${PROJECT_BINARY_DIR}/DocProjectFileAux)
set(DocumentationProjectFileInputDir ${PROJECT_SOURCE_DIR}/Documentation/ProjectFile)

# remove old output
if (IS_DIRECTORY ${DocumentationProjectFileBuildDir})
	file(REMOVE_RECURSE ${DocumentationProjectFileBuildDir})
endif()

file(GLOB_RECURSE input_paths ${DocumentationProjectFileInputDir}/_*)

foreach(p ${input_paths})
	message("in path ${p}")
	documentationProjectFilePutIntoPlace(${p})
endforeach()
