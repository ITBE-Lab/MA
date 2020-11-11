@ECHO OFF
rem See: https://stackoverflow.com/questions/4580101/python-add-pythonpath-during-command-line-module-run
rem %1: Python-Executable
rem %2: PYTHONPATH
rem %3: Filename of Python-Script

setlocal
set PYTHONPATH=.
rem :%PYTHONPATH%
echo "Python Interpreter: %1"
echo "PYTHONPATH %PYTHONPATH%"
%1 %3
endlocal