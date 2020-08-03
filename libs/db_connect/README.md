## Notes
### MySQL
So far, there comes no FindMYSQL with the cmake:  
Got some from the WEB: https://gitlab.kitware.com/vtk/vtk/blob/master/CMake/FindMySQL.cmake  
In this file you have to insertion your release number to _MySQL_versions or 

  
| Filename | Description |  
|---|---|
| db_sql.h | Primary database interface for MA for now  | 

# Optimizing Server settings:

innodb_buffer_pool_size = 4G
concurrent_insert = ALWAYS
datadir = /mnt/ssd1/mysql                       # original: /var/lib/mysql
log-error = /mnt/ssd1/mysql-error-log/error.log
tmpdir = /mnt/ssd1/tmp
