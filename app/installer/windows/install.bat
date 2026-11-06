@echo off
rem Installer of the eHabitat+ app for Windows (see the page Desktop app of the wiki of eHabpy). It is opened with a
rem double click from the folder app\installer\windows of a copy of the repository (Windows may warn: "More info" >
rem "Run anyway"), or downloaded and opened with a command in the Command Prompt (what curl downloads has no mark of
rem the Internet, so SmartScreen does not stop it):
rem
rem   curl -fsSLo "%TEMP%\install-ehabitat.bat" https://raw.githubusercontent.com/javimarlop/eHabpy/master/app/installer/windows/install.bat && "%TEMP%\install-ehabitat.bat"
rem
rem It installs in %USERPROFILE%\eHabitat, without administrator rights:
rem   - micromamba (the package installer of conda-forge), in %USERPROFILE%\eHabitat\micromamba;
rem   - the environment of eHabitat+ (Python, GDAL, GRASS GIS, R and its packages, Earth Engine, MonitorEO and
rem     pywebview, from the environment.yml of the repository), in %USERPROFILE%\eHabitat\env;
rem   - the code of eHabitat+ in %USERPROFILE%\eHabitat\eHabpy, only if the installer is not opened from a copy of
rem     the repository (if it is, the app uses that copy as its working folder).
rem It puts on the Desktop the shortcut "eHabitat+", which opens the app in its window with pythonw.exe, without a
rem console (its messages go to %USERPROFILE%\eHabitat\ehabitat.log), and opens it. If the environment already exists
rem and is complete, it only creates the shortcut; to install it again, delete the folder env or set EHAB_REINSTALL=1.
rem Other variables: EHAB_HOME (installation folder), EHAB_BRANCH (branch of GitHub of the code it downloads) and
rem EHAB_NO_OPEN=1 (do not open the app at the end). No accents: cmd reads this file with the code page of the console.
setlocal EnableExtensions
title Installer of eHabitat+
set "REPO=javimarlop/eHabpy"
if not defined EHAB_BRANCH set "EHAB_BRANCH=master"
if defined EHAB_HOME (set "DIR=%EHAB_HOME%") else (set "DIR=%USERPROFILE%\eHabitat")
set "ENV=%DIR%\env"
set "MAMBA_DIR=%DIR%\micromamba"
set "MAMBA=%MAMBA_DIR%\micromamba.exe"
set "PYTHONNOUSERSITE=1"

echo Installer of eHabitat+
echo ======================
echo.
echo It installs in %DIR%, without administrator rights, what eHabitat+ needs:
echo Python, GDAL, GRASS GIS, R, Google Earth Engine and MonitorEO (about 0.7 GB to download and 3 GB on disk,
echo only the first time).
echo.
if not exist "%DIR%" mkdir "%DIR%" || goto :fail

rem Code: the copy of the repository that contains this installer or, if none, that of %DIR%\eHabpy
set "HERE=%~dp0"
if exist "%HERE%..\..\ehabitat.py" if exist "%HERE%..\..\..\environment.yml" (
  for %%I in ("%HERE%..\..\..") do set "CODE=%%~fI"
  goto :code_ok
)
set "CODE=%DIR%\eHabpy"
if exist "%CODE%\app\ehabitat.py" (
  echo Code of eHabitat+: %CODE% ^(already downloaded; left as it is^)
  goto :code_ok
)
echo Downloading the code of eHabitat+ (github.com/%REPO%, branch %EHAB_BRANCH%)...
set "TMPD=%DIR%\download"
if exist "%TMPD%" rmdir /s /q "%TMPD%"
mkdir "%TMPD%" || goto :fail
curl -fL --progress-bar -o "%TMPD%\eHabpy.zip" "https://github.com/%REPO%/archive/refs/heads/%EHAB_BRANCH%.zip" || goto :fail
tar -xf "%TMPD%\eHabpy.zip" -C "%TMPD%" || goto :fail
for /d %%D in ("%TMPD%\*") do set "SRC=%%~fD"
if not exist "%SRC%\app\ehabitat.py" (
  echo The branch %EHAB_BRANCH% of %REPO% has no app ^(app\ehabitat.py^).
  goto :fail
)
move "%SRC%" "%CODE%" >nul || goto :fail
rmdir /s /q "%TMPD%"
:code_ok
echo Code of eHabitat+: %CODE%

if "%EHAB_REINSTALL%"=="1" goto :install
call :env_ok && (
  echo The environment is already installed in %ENV%.
  goto :launcher
)

:install
if not exist "%MAMBA_DIR%" mkdir "%MAMBA_DIR%"
if exist "%MAMBA%" (
  echo 1/3 micromamba is already in %MAMBA_DIR%.
  goto :mamba_ok
)
echo 1/3 Downloading micromamba (the package installer of conda-forge)...
set "URL=https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-win-64"
curl -fL --progress-bar -o "%MAMBA%.part" "%URL%" || goto :fail
curl -fsL -o "%MAMBA%.sha256" "%URL%.sha256" || goto :fail
powershell -NoProfile -Command "if ((Get-FileHash -Algorithm SHA256 -LiteralPath ($env:MAMBA + '.part')).Hash -ne (Get-Content -Raw -LiteralPath ($env:MAMBA + '.sha256')).Trim()) { exit 1 }" || goto :damaged
move /y "%MAMBA%.part" "%MAMBA%" >nul || goto :fail
:mamba_ok

echo 2/3 List of packages: %CODE%\environment.yml
echo 3/3 Installing Python, GDAL, GRASS GIS, R, Earth Engine and MonitorEO (a few minutes)...
if exist "%ENV%" rmdir /s /q "%ENV%"
set "MAMBA_ROOT_PREFIX=%MAMBA_DIR%"
"%MAMBA%" create -y -p "%ENV%" -f "%CODE%\environment.yml" || goto :fail
"%MAMBA%" clean -y --all >nul 2>&1
call :env_ok || (
  echo The environment was installed, but something is missing ^(Python, R, GRASS GIS or pywebview^).
  goto :fail
)

:launcher
rem Shortcut of the Desktop: pythonw.exe (no console) with the app; the app activates the environment by itself
set "ICON=%CODE%\app\ehabitat_app\ui\icon.ico"
if not exist "%ICON%" set "ICON=%ENV%\python.exe"
powershell -NoProfile -Command "$s = (New-Object -ComObject WScript.Shell).CreateShortcut((Join-Path ([Environment]::GetFolderPath('Desktop')) 'eHabitat+.lnk')); $s.TargetPath = Join-Path $env:ENV 'pythonw.exe'; $s.Arguments = '-s \"' + (Join-Path $env:CODE 'app\ehabitat.py') + '\"'; $s.WorkingDirectory = $env:CODE; $s.IconLocation = $env:ICON + ',0'; $s.Description = 'eHabitat+'; $s.Save()" || goto :fail
rem Launcher with a console, to see the messages if the window does not open
> "%DIR%\eHabitat+ (with console).bat" (
  echo @echo off
  echo set "PYTHONNOUSERSITE=1"
  echo set "EHAB_MAMBA=%MAMBA%"
  echo cd /d "%CODE%"
  echo "%ENV%\python.exe" -s "%CODE%\app\ehabitat.py" %%*
  echo pause
)

echo.
echo Done. To open eHabitat+: double click "eHabitat+" on the Desktop. You can close this window now.
echo If the window of the app does not open, use "%DIR%\eHabitat+ (with console).bat" to see the messages,
echo or open the app in the web browser with: "%DIR%\eHabitat+ (with console).bat" --browser
if not "%EHAB_NO_OPEN%"=="1" (
  set "EHAB_MAMBA=%MAMBA%"
  start "" /d "%CODE%" "%ENV%\pythonw.exe" -s "%CODE%\app\ehabitat.py"
)
exit /b 0

:env_ok
if not exist "%ENV%\python.exe" exit /b 1
if not exist "%ENV%\pythonw.exe" exit /b 1
if not exist "%ENV%\Library\bin\grass.bat" exit /b 1
if not exist "%ENV%\Scripts\Rscript.exe" if not exist "%ENV%\lib\R\bin\Rscript.exe" exit /b 1
"%ENV%\python.exe" -s -c "import webview, osgeo.gdal, yaml, geopandas, esda, ee" >nul 2>&1 || exit /b 1
exit /b 0

:damaged
echo.
echo ERROR: the download of micromamba is damaged.
:fail
echo.
echo ERROR: the installation did not finish. Check the internet connection and open the installer again.
echo The messages above say what failed.
pause
exit /b 1
