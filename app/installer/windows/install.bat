@echo off
rem Installer of the eHabitat+ app for Windows (see the page Desktop app of the wiki of eHabpy). It is opened with a
rem double click from the folder app\installer\windows of a copy of the repository (Windows may warn: "More info" >
rem "Run anyway"), or downloaded and opened with a command in the Command Prompt (what curl downloads has no mark of
rem the Internet, so SmartScreen does not stop it):
rem
rem   curl -fsSLo "%TEMP%\install-ehabitat.bat" https://raw.githubusercontent.com/javimarlop/eHabpy/master/app/installer/windows/install.bat && "%TEMP%\install-ehabitat.bat"
rem
rem It installs in %USERPROFILE%\eHabitat, without administrator rights:
rem   - micromamba (the package installer of conda-forge), in %USERPROFILE%\eHabitat\micromamba, with the DLLs of the
rem     runtime of Visual C++ next to it when this computer does not have that runtime (a new Windows);
rem   - the environment of eHabitat+ (Python, GDAL, GRASS GIS, R and its packages, Earth Engine, MonitorEO and
rem     pywebview, from the environment.yml of the repository), in %USERPROFILE%\eHabitat\env;
rem   - the code of eHabitat+ in %USERPROFILE%\eHabitat\eHabpy, only if the installer is not opened from a copy of
rem     the repository (if it is, the app uses that copy as its working folder).
rem It puts on the Desktop the shortcut "eHabitat+", which opens the app in its window with pythonw.exe, without a
rem console (its messages go to %USERPROFILE%\eHabitat\ehabitat.log), and opens it. If the environment already exists
rem and is complete, it only checks the R packages that come from CRAN and creates the shortcut; to install it again,
rem delete the folder env or set EHAB_REINSTALL=1.
rem The app keeps itself updated from GitHub every time it opens (see Updates in that page of the wiki).
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
  goto :code_shown
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
:code_shown

rem Packages: those of environment.yml but the ones marked [win: CRAN name] there, which conda-forge does not have for
rem Windows; those (R packages) are installed from CRAN after the environment
set "ENVFILE=%DIR%\environment-windows.yml"
set "CRANLIST=%DIR%\cran-packages.txt"
powershell -NoProfile -Command "$cran = @(); $out = foreach ($l in [IO.File]::ReadAllLines((Join-Path $env:CODE 'environment.yml'))) { if ($l -match '\[win: CRAN ([A-Za-z0-9.]+)\]') { $cran += $Matches[1] } else { $l } }; $u = New-Object Text.UTF8Encoding($false); [IO.File]::WriteAllLines($env:ENVFILE, [string[]]$out, $u); [IO.File]::WriteAllLines($env:CRANLIST, [string[]]$cran, $u)" || goto :fail

if "%EHAB_REINSTALL%"=="1" goto :install
call :openmp
call :env_ok && (
  echo The environment is already installed in %ENV%.
  goto :cran_step
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
rem micromamba needs the runtime of Visual C++ (MSVCP140.dll, VCRUNTIME140.dll, ...), which Windows does not have until
rem a program installs it, and whose installer needs administrator rights: if micromamba does not open, the DLLs of the
rem package vc14_runtime of conda-forge go next to it, where Windows looks for them first
"%MAMBA%" --version >nul 2>&1 || call :vcruntime || goto :fail

echo 2/3 List of packages: %CODE%\environment.yml
echo 3/3 Installing Python, GDAL, GRASS GIS, R, Earth Engine and MonitorEO (a few minutes)...
if exist "%ENV%" rmdir /s /q "%ENV%"
set "MAMBA_ROOT_PREFIX=%MAMBA_DIR%"
rem Cache of the packages: micromamba keeps them in folders like pkgs\https\conda.anaconda.org\conda-forge\noarch\<package>
rem and, from %DIR%, the C++ headers of the compilers that R needs go over the limit of 260 characters of the paths of
rem Windows: the cache goes to a new short folder of the disk, deleted at the end. Only this user (and SYSTEM and the
rem administrators, as in the profile) can open it: a folder made in C:\ can be changed by every user of the computer,
rem and micromamba links the files of the cache into the environment, where they keep the permissions of the cache
set "CONDA_PKGS_DIRS=%SystemDrive%\ehp-%RANDOM%%RANDOM%"
mkdir "%CONDA_PKGS_DIRS%" 2>nul || set "CONDA_PKGS_DIRS=%DIR%\pkgs"
if /i not "%CONDA_PKGS_DIRS%"=="%DIR%\pkgs" icacls "%CONDA_PKGS_DIRS%" /inheritance:r /grant:r "%USERDOMAIN%\%USERNAME%:(OI)(CI)F" "*S-1-5-18:(OI)(CI)F" "*S-1-5-32-544:(OI)(CI)F" >nul 2>&1 || (
  rmdir "%CONDA_PKGS_DIRS%" 2>nul
  set "CONDA_PKGS_DIRS=%DIR%\pkgs"
)
"%MAMBA%" create -y -p "%ENV%" -f "%ENVFILE%" || (
  rmdir /s /q "%CONDA_PKGS_DIRS%" 2>nul
  goto :fail
)
"%MAMBA%" clean -y --all >nul 2>&1
rmdir /s /q "%CONDA_PKGS_DIRS%" 2>nul
call :openmp
call :env_ok || (
  echo The environment was installed, but something is missing ^(Python, R, GRASS GIS or pywebview^).
  goto :fail
)
:cran_step
call :cran || goto :fail
call :grass_db || (
  echo.
  echo WARNING: GRASS GIS cannot open its databases ^(the attributes of the vector maps^), a fault of its build 8.5.0
  echo          of conda-forge for Windows: the HFT model cannot run here until conda-forge publishes a fixed build
  echo          ^(then delete the folder env and open this installer again^). Until then use eHabitat+ in WSL2,
  echo          with its installer for Ubuntu: see Windows in the page Desktop app of the wiki.
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
  set "CONDA_PKGS_DIRS="
  start "" /d "%CODE%" "%ENV%\pythonw.exe" -s "%CODE%\app\ehabitat.py"
)
exit /b 0

:env_ok
if not exist "%ENV%\python.exe" exit /b 1
if not exist "%ENV%\pythonw.exe" exit /b 1
if not exist "%ENV%\Library\bin\grass.bat" exit /b 1
if not exist "%ENV%\Scripts\Rscript.exe" if not exist "%ENV%\lib\R\bin\Rscript.exe" exit /b 1
"%ENV%\python.exe" -s -c "import webview, osgeo.gdal, yaml, geopandas, esda, ee" >nul 2>&1 || exit /b 1
rem A module of GRASS GIS opens, in a temporary project (grass.bat changes variables: only here, with setlocal)
setlocal
set "PATH=%ENV%;%ENV%\Library\mingw-w64\bin;%ENV%\Library\usr\bin;%ENV%\Library\bin;%ENV%\Scripts;%ENV%\bin;%PATH%"
call "%ENV%\Library\bin\grass.bat" --tmp-project EPSG:4326 --exec g.proj -p >nul 2>&1 || exit /b 1
exit /b 0

:grass_db
rem The databases of GRASS GIS (the attributes of the vector maps, which the HFT model needs) open: the build 8.5.0 of
rem GRASS of conda-forge for Windows fails there (its drivers call the empty functions of grass_dbstubs.dll)
setlocal
set "PATH=%ENV%;%ENV%\Library\mingw-w64\bin;%ENV%\Library\usr\bin;%ENV%\Library\bin;%ENV%\Scripts;%ENV%\bin;%PATH%"
call "%ENV%\Library\bin\grass.bat" --tmp-project EPSG:4326 --exec db.execute "sql=create table t (a integer)" >nul 2>&1 || exit /b 1
exit /b 0

:openmp
rem GRASS GIS 8.5 of conda-forge for Windows is built with the OpenMP of LLVM of Visual Studio, libomp140.x86_64.dll,
rem which Microsoft does not redistribute and Windows does not have: without it no module of GRASS opens. The package
rem llvm-openmp has that runtime (libomp.dll) and libiomp5md.dll, which only forwards its functions to libomp.dll: a
rem copy of it with the name that GRASS asks for gives GRASS that runtime, the same one of the rest of the environment
if exist "%ENV%\Library\bin\libomp140.x86_64.dll" exit /b 0
if exist "%ENV%\Library\bin\libiomp5md.dll" copy /y "%ENV%\Library\bin\libiomp5md.dll" "%ENV%\Library\bin\libomp140.x86_64.dll" >nul
exit /b 0

:cran
rem The R packages of CRANLIST that R cannot load, from CRAN (binary packages of Windows)
if not exist "%CRANLIST%" exit /b 0
for %%F in ("%CRANLIST%") do if %%~zF==0 exit /b 0
setlocal
set "PATH=%ENV%;%ENV%\Library\mingw-w64\bin;%ENV%\Library\usr\bin;%ENV%\Library\bin;%ENV%\Scripts;%ENV%\bin;%PATH%"
set "RSCRIPT=%ENV%\Scripts\Rscript.exe"
if not exist "%RSCRIPT%" set "RSCRIPT=%ENV%\lib\R\bin\Rscript.exe"
"%RSCRIPT%" -e "p <- readLines(Sys.getenv('CRANLIST')); p <- p[nzchar(p)]; m <- p[!sapply(p, requireNamespace, quietly = TRUE)]; if (length(m)) { message('   From CRAN (conda-forge does not have them for Windows): ', paste(m, collapse = ', ')); install.packages(m, repos = 'https://cloud.r-project.org', type = 'binary', quiet = TRUE); m <- m[!sapply(m, requireNamespace, quietly = TRUE)]; if (length(m)) { message('R packages not installed: ', paste(m, collapse = ', ')); quit(status = 1) } }" || exit /b 1
exit /b 0

:vcruntime
echo    micromamba needs the runtime of Visual C++, which is not on this computer: downloading it from conda-forge...
powershell -NoProfile -Command "$ErrorActionPreference = 'Stop'; $ProgressPreference = 'SilentlyContinue'; [Net.ServicePointManager]::SecurityProtocol = [Net.ServicePointManager]::SecurityProtocol -bor 3072; $all = Invoke-RestMethod 'https://api.anaconda.org/package/conda-forge/vc14_runtime/files'; $f = $all | Where-Object { $_.attrs.subdir -eq 'win-64' -and $_.labels -contains 'main' -and $_.basename -like '*.conda' } | Sort-Object { [version]$_.version }, { [int]$_.attrs.build_number } | Select-Object -Last 1; $t = Join-Path $env:MAMBA_DIR 'vc14_runtime'; if (Test-Path $t) { Remove-Item -Recurse -Force $t }; New-Item -ItemType Directory $t | Out-Null; $p = Join-Path $t 'package.conda'; Invoke-WebRequest -UseBasicParsing ('https:' + $f.download_url) -OutFile $p; if ((Get-FileHash -Algorithm SHA256 -LiteralPath $p).Hash -ne $f.sha256) { throw 'the download is damaged' }; $tar = Join-Path $env:SystemRoot 'System32\tar.exe'; & $tar -xf $p -C $t; if ($LASTEXITCODE) { throw 'tar could not open the package' }; $z = Get-ChildItem $t -Filter 'pkg-*.tar.zst' | Select-Object -First 1; & $tar -xf $z.FullName -C $t; if ($LASTEXITCODE) { throw 'tar could not open the package' }; Copy-Item (Join-Path $t '*.dll') $env:MAMBA_DIR -Force; Remove-Item -Recurse -Force $t; Write-Host ('   ' + $f.basename)" || goto :vcruntime_fail
"%MAMBA%" --version >nul 2>&1 && exit /b 0
echo micromamba does not open, not even with the runtime of Visual C++ next to it.
:vcruntime_fail
rem e.g. Windows 10, whose tar does not open the .tar.zst inside the packages of conda-forge
echo Install the runtime of Visual C++ ^(x64; it needs administrator rights^) from
echo https://aka.ms/vs/17/release/vc_redist.x64.exe and open this installer again.
exit /b 1

:damaged
echo.
echo ERROR: the download of micromamba is damaged.
:fail
echo.
echo ERROR: the installation did not finish. Check the internet connection and open the installer again.
echo The messages above say what failed.
pause
exit /b 1
