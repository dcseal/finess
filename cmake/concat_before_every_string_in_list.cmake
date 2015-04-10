function(concat_before_every_string_in_list head output)
    foreach(p ${ARGN})
        list(APPEND new "${head}${p}")
    endforeach()
    set(${output} ${new} PARENT_SCOPE)
endfunction()
