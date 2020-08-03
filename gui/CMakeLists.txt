# CMakeList.txt : CMake for wxWidgets based GUI of MA
# Builds the MA GUI-application
# Requires a wxWdigest build with static libraries (-DwxBUILD_SHARED=OFF)

cmake_minimum_required (VERSION 3.8)

project( "MAgui" VERSION 0.1 DESCRIPTION "Graphical User Interface for MA" )

# Find wxWidgets includes a libraries
# Hint: 'adv' integrates the wxAdvanced library 
#		 (that advanced library is required because of wxWizard) 
find_package( wxWidgets COMPONENTS core base adv )

if(wxWidgets_FOUND)
	message( "Include GUI app..." )
	# Define a variable containing a list of source files for the project
 
	set( SRC_FILES src/maGUI.cpp res/maGUI.rc )

	if(MSVC)
		add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	endif()
	
	# Path to wxWidges includes
	include( ${wxWidgets_USE_FILE} )

	# Use static binding (wxWidgets must support this as well!)
	# (No wxWidget DLLs required during runtime.)
	set(wxWidgets_USE_STATIC ON)

	# Build the MA executable itself
	add_executable ( maGUI WIN32 MACOSX_BUNDLE ${SRC_FILES} )

	target_link_libraries( maGUI libMA ${wxWidgets_LIBRARIES} )
    target_compile_definitions( maGUI PUBLIC DLL_PORT_MA=dllimport )

	# Local includes 
	target_include_directories( maGUI PUBLIC inc )
	
	# build maGUI.exe in the main folder
	set_target_properties(maGUI PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )
endif()
