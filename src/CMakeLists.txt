cmake_minimum_required(VERSION 3.16.0)

set(PROJECT_NAME isChIP)
project(${PROJECT_NAME})

# add_compile_definitions(
	# _ISCHIP
	# _FEATURES
	# _FILE_WRITE
	# _MULTITHREAD
	# )
	
# ZLIB
if(MSVC)
	set(ZLIB_DIR /ProgramData/zlib1.3)
	set(ZLIB_LIBRARY_DEBUG ${ZLIB_DIR}/Debug/zlibstaticd.lib)
	set(ZLIB_LIBRARY_RELEASE ${ZLIB_DIR}/Release/zlibstatic.lib)
	set(ZLIB_INCLUDE_DIRS ${ZLIB_DIR})
endif()

find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
	# add_compile_definitions(_ZLIB)
	include_directories( ${ZLIB_INCLUDE_DIRS} )
else()
	message(WARNING ">> ZLIB not found!")
endif()

# Source groups
set(BASIC_DIR .)
# set(BASIC_DIR ../../bioinfoBasic)

# include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/${BASIC_DIR} )

set(Header_Files
	${BASIC_DIR}/common.h
	${BASIC_DIR}/Data.h
	${BASIC_DIR}/DataInFile.h
	${BASIC_DIR}/TxtFile.h
	isChIP.h
	Imitator.h
	DataOutFile.h
	effPartition.h
	RandomGen.h
)
source_group("Header Files" FILES ${Header_Files})

set(Source_Files
	${BASIC_DIR}/common.cpp
	${BASIC_DIR}/Data.cpp
	${BASIC_DIR}/DataInFile.cpp
	${BASIC_DIR}/TxtFile.cpp
	isChIP.cpp
	Imitator.cpp
	DataOutFile.cpp
	effPartition.cpp
	RandomGen.cpp
)

source_group("Source Files" FILES ${Source_Files})

set(ALL_FILES
	${Header_Files}
	${Source_Files}
)

# Target
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY build)
add_executable(${PROJECT_NAME} ${ALL_FILES})

if(ZLIB_FOUND)
	target_link_libraries(${PROJECT_NAME} PRIVATE ZLIB::ZLIB)
endif()

set(ROOT_NAMESPACE ${PROJECT_NAME})

set_target_properties(${PROJECT_NAME} PROPERTIES
	VS_GLOBAL_KEYWORD "Win32Proj"
)

if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
	set_target_properties(${PROJECT_NAME} PROPERTIES
		INTERPROCEDURAL_OPTIMIZATION_RELEASE "TRUE"
	)
elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
	set_target_properties(${PROJECT_NAME} PROPERTIES
		INTERPROCEDURAL_OPTIMIZATION_DEBUG   "FALSE"
		INTERPROCEDURAL_OPTIMIZATION_RELEASE "TRUE"
	)
endif()

# Compile definitions
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
	target_compile_definitions(${PROJECT_NAME} PRIVATE
		"$<$<CONFIG:Debug>:"
			"_DEBUG"
		">"
		"$<$<CONFIG:Release>:"
			"NDEBUG"
		">"
		"_CONSOLE;"
		"UNICODE;"
		"_UNICODE"
		"_CRT_SECURE_NO_WARNINGS"
	)
elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
	target_compile_definitions(${PROJECT_NAME} PRIVATE
		"$<$<CONFIG:Debug>:"
			"_DEBUG"
		">"
		"$<$<CONFIG:Release>:"
			"NDEBUG"
		">"
		"WIN32;"
		"_CONSOLE;"
		"_CRT_SECURE_NO_WARNINGS"
	)
endif()

# Compile and link options
if(MSVC)
	if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
		target_compile_options(${PROJECT_NAME} PRIVATE
			$<$<CONFIG:Release>:
				/Oi;
				/Gy
			>
			/permissive-;
			/sdl;
			/W3;
			${DEFAULT_CXX_DEBUG_INFORMATION_FORMAT};
			${DEFAULT_CXX_EXCEPTION_HANDLING}
		)
	elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
		target_compile_options(${PROJECT_NAME} PRIVATE
			$<$<CONFIG:Debug>:
				/JMC;
				/Od;
				/Oi-
			>
			$<$<CONFIG:Release>:
				/Oi;
				/Gy
			>
			/permissive-;
			/std:c++17;
			/sdl;
			/W3;
			${DEFAULT_CXX_DEBUG_INFORMATION_FORMAT};
			${DEFAULT_CXX_EXCEPTION_HANDLING}
		)
	endif()
	if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
		target_link_options(${PROJECT_NAME} PRIVATE
			$<$<CONFIG:Debug>:
				/INCREMENTAL
			>
			$<$<CONFIG:Release>:
				/OPT:REF;
				/OPT:ICF;
				/INCREMENTAL:NO
			>
			/DEBUG;
			/SUBSYSTEM:CONSOLE
		)
	elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
		target_link_options(${PROJECT_NAME} PRIVATE
			$<$<CONFIG:Release>:
				/OPT:ICF
			>
			/DEBUG;
			/SUBSYSTEM:CONSOLE;
			/OPT:REF;
			/INCREMENTAL:NO
		)
	endif()
endif()
