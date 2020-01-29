# - Try to find MySQL.
# Once done this will define:
# MYSQL_FOUND			- If false, do not try to use MySQL.
# MYSQL_INCLUDE_DIRS	- Where to find mysql.h, etc.
# MYSQL_LIBRARIES		- The libraries to link against.
# MYSQL_VERSION_STRING	- Version in a string of MySQL.
#
# Created by RenatoUtsch based on eAthena implementation.
#
# Please note that this module only supports Windows and Linux officially, but
# should work on all UNIX-like operational systems too.
#

#=============================================================================
# Copyright 2012 RenatoUtsch
# Copyright 2018 Nicolas Mora <mail@babelouest.org>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)
# "$ENV{PROGRAMFILES(x86)}/MySQL/*/include"
# "$ENV{PROGRAMFILES(x86)}/MySQL/*/lib"
if( WIN32 )
	set(_MySQL_mariadb_versions 10.2 10.3)
	  set(_MySQL_versions 5.0 5.7)
	  set(_MySQL_paths)
	  foreach (_MySQL_version IN LISTS _MySQL_mariadb_versions)
	    list(APPEND _MySQL_paths
	      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MariaDB ${_MySQL_version};INSTALLDIR]"
	      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MariaDB ${_MySQL_version} (x64);INSTALLDIR]")
	  endforeach ()
	  foreach (_MySQL_version IN LISTS _MySQL_versions)
	    list(APPEND _MySQL_paths
	      "C:/Program Files/MySQL/MySQL Server ${_MySQL_version}/lib/opt"
	      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MySQL AB\\MySQL Server ${_MySQL_version};Location]"
	      "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Wow6432Node\\MySQL AB\\MySQL Server ${_MySQL_version};Location]")
	  endforeach ()
	  unset(_MySQL_version)
	  unset(_MySQL_versions)
	  unset(_MySQL_mariadb_versions)
	
	  find_path(MySQL_INCLUDE_DIR
	    NAMES mysql.h
	    PATHS
	      "C:/Program Files/MySQL/include"
	      "C:/MySQL/include"
	      ${_MySQL_paths}
	    PATH_SUFFIXES include include/mysql
	    DOC "Location of mysql.h")
	  mark_as_advanced(MySQL_INCLUDE_DIR)
	  find_library(MySQL_LIBRARY
	    NAMES libmariadb mysql libmysql mysqlclient
	    PATHS
	      "C:/Program Files/MySQL/lib"
	      "C:/MySQL/lib/debug"
	      ${_MySQL_paths}
	    PATH_SUFFIXES lib lib/opt
	    DOC "Location of the mysql library")
	  mark_as_advanced(MySQL_LIBRARY) 
	
	  include(FindPackageHandleStandardArgs)
	  find_package_handle_standard_args(MySQL
	    REQUIRED_VARS MySQL_INCLUDE_DIR MySQL_LIBRARY)
	
	  if (MySQL_FOUND)
	    set(MYSQL_INCLUDE_DIR "${MySQL_INCLUDE_DIR}")
	    set(MYSQL_LIBRARY "${MySQL_LIBRARY}")
	    if (NOT TARGET MySQL::MySQL)
	      add_library(MySQL::MySQL UNKNOWN IMPORTED)
	      set_target_properties(MySQL::MySQL PROPERTIES
	      IMPORTED_LOCATION "${MySQL_LIBRARY}"
	        INTERFACE_INCLUDE_DIRECTORIES "${MySQL_INCLUDE_DIR}")
	    endif ()
	  endif ()
		# find_path( MYSQL_INCLUDE_DIR
		# 	NAMES "mysql.h"
		# 	PATHS "C:/Program Files/MySQL/*/include"
		# 		  "$ENV{SYSTEMDRIVE}/MySQL/*/include" )
		# 
		# find_library( MYSQL_LIBRARY
		# 	NAMES "mysqlclient" "mysqlclient_r"
		# 	PATHS "C:/Program Files/MySQL/*/lib"
		# 		  "$ENV{SYSTEMDRIVE}/MySQL/*/lib" )
else()
	find_path( MYSQL_INCLUDE_DIR
		NAMES "mysql.h"
		PATHS "/usr/include/mysql"
			  "/usr/local/include/mysql"
			  "/usr/mysql/include/mysql" )
	
	find_library( MYSQL_LIBRARY
		NAMES "mysqlclient" "mysqlclient_r" "mariadbclient" "mariadbclient_r" 
		PATHS "/usr/lib/x86_64-linux-gnu"
			  "/lib/mysql"
			  "/lib/arm-linux-gnueabihf"
			  "/lib64/mysql"
			  "/usr/lib/mysql"
			  "/usr/lib64/mysql"
			  "/usr/lib/arm-linux-gnueabihf"
			  "/usr/local/lib/mysql"
			  "/usr/local/lib64/mysql"
			  "/usr/local/lib/arm-linux-gnueabihf"
			  "/usr/mysql/lib/mysql"
			  "/usr/mysql/lib64/mysql"
			  ) 
endif()

set(MYSQL_VERSION_STRING "0.0.0")
if( MYSQL_INCLUDE_DIR AND EXISTS "${MYSQL_INCLUDE_DIR}/mysql_version.h" )
    set(regex_version "^#define[ \t]+MYSQL_SERVER_VERSION[ \t]+\"([^\"]+)\".*")
    file(STRINGS "${MYSQL_INCLUDE_DIR}/mysql_version.h" mysql_version REGEX "${regex_version}")
    string(REGEX REPLACE "${regex_version}" "\\1" MYSQL_VERSION_STRING "${mysql_version}")
    unset(regex_version)
    unset(mysql_version)
endif()

# handle the QUIETLY and REQUIRED arguments and set MYSQL_FOUND to TRUE if
# all listed variables are TRUE
include( FindPackageHandleStandardArgs )
find_package_handle_standard_args( MYSQL REQUIRED_VARS	
  MYSQL_LIBRARY MYSQL_INCLUDE_DIR
	VERSION_VAR		MYSQL_VERSION_STRING )

set( MYSQL_INCLUDE_DIRS ${MYSQL_INCLUDE_DIR} )
set( MYSQL_LIBRARIES ${MYSQL_LIBRARY} )

message( "FindMySQL - Include:" ${MYSQL_INCLUDE_DIRS} " Library: " ${MYSQL_LIBRARIES} )

mark_as_advanced( MYSQL_INCLUDE_DIR MYSQL_LIBRARY )
