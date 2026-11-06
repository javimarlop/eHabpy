#!/bin/bash
# Installer of the eHabitat+ app for macOS (see the page Desktop app of the wiki of eHabpy). It is opened with a
# double click from the folder app/installer/macos of a copy of the repository (the first time macOS may ask for
# permission: right click > Open), or pasted in the Terminal:
#
#   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/javimarlop/eHabpy/master/app/installer/macos/install.command)"
#
# It installs in ~/eHabitat, without administrator rights:
#   - micromamba (the package installer of conda-forge), in ~/eHabitat/micromamba;
#   - the environment of eHabitat+ (Python, GDAL, GRASS GIS, R and its packages, Earth Engine, MonitorEO and
#     pywebview, from the environment.yml of the repository), in ~/eHabitat/env;
#   - the code of eHabitat+ in ~/eHabitat/eHabpy, only if the installer is not opened from a copy of the repository
#     (if it is, the app uses that copy as its working folder: inVars, segm_optim/pa, results...).
# It puts on the Desktop the app "eHabitat+", which opens the app in its window, without the Terminal (its messages
# go to ~/eHabitat/ehabitat.log), and opens it. If the environment already exists and is complete, it only creates
# the app; to install it again, delete ~/eHabitat/env or open the installer with EHAB_REINSTALL=1.
# Other variables: EHAB_HOME (installation folder), EHAB_DESKTOP (folder of the shortcut), EHAB_BRANCH (branch of
# GitHub of the code it downloads) and EHAB_NO_OPEN=1 (do not open the app at the end).
# Everything is inside main, which is called at the end: a broken download runs nothing.

main() {
	set -euo pipefail
	REPO=javimarlop/eHabpy
	BRANCH="${EHAB_BRANCH:-master}"
	DIR="${EHAB_HOME:-$HOME/eHabitat}"
	ENV="$DIR/env"
	MAMBA_DIR="$DIR/micromamba"
	MAMBA="$MAMBA_DIR/micromamba"
	DESKTOP="${EHAB_DESKTOP:-$HOME/Desktop}"
	APP="$DESKTOP/eHabitat+.app"

	fail() {
		echo
		echo "ERROR: $1"
		echo "Check the internet connection and open the installer again. The messages above say what failed."
		read -r -p "Press Enter to close. " _ < /dev/tty || true
		exit 1
	}
	trap 'fail "the installation stopped at line $LINENO."' ERR

	echo "Installer of eHabitat+"
	echo "======================"
	echo
	echo "It installs in $DIR, without administrator rights, what eHabitat+ needs: Python, GDAL, GRASS GIS, R,"
	echo "Google Earth Engine and MonitorEO (about 0.6 GB to download and 2.5 GB on disk, only the first time)."
	echo

	case "$(uname -m)" in
		arm64) PLATFORM=osx-arm64 ;;
		x86_64) PLATFORM=osx-64 ;;
		*) fail "this Mac is neither Apple Silicon nor Intel ($(uname -m))." ;;
	esac
	mkdir -p "$DIR"

	# Code: the copy of the repository that contains this installer or, if none, that of ~/eHabitat/eHabpy
	HERE=""
	if [ -n "${BASH_SOURCE[0]:-}" ] && [ -f "${BASH_SOURCE[0]}" ]; then
		HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
	fi
	if [ -n "$HERE" ] && [ -f "$HERE/../../ehabitat.py" ] && [ -f "$HERE/../../../environment.yml" ]; then
		CODE="$(cd "$HERE/../../.." && pwd)"
		echo "Code of eHabitat+: $CODE"
	else
		CODE="$DIR/eHabpy"
		if [ -f "$CODE/app/ehabitat.py" ]; then
			echo "Code of eHabitat+: $CODE (already downloaded; left as it is)"
		else
			echo "Downloading the code of eHabitat+ (github.com/$REPO, branch $BRANCH)..."
			TMP="$DIR/download"
			rm -rf "$TMP" && mkdir -p "$TMP"
			curl -fL --progress-bar -o "$TMP/eHabpy.zip" "https://github.com/$REPO/archive/refs/heads/$BRANCH.zip" \
				|| fail "the code of eHabitat+ could not be downloaded."
			(cd "$TMP" && unzip -q eHabpy.zip) || fail "the downloaded code is damaged."
			SRC="$(find "$TMP" -maxdepth 1 -mindepth 1 -type d | head -n 1)"
			[ -f "$SRC/app/ehabitat.py" ] || fail "the branch $BRANCH of $REPO has no app (app/ehabitat.py)."
			mv "$SRC" "$CODE"
			rm -rf "$TMP"
		fi
	fi

	env_ok() {
		[ -x "$ENV/bin/python" ] && [ -x "$ENV/bin/Rscript" ] && [ -x "$ENV/bin/grass" ] \
			&& PYTHONNOUSERSITE=1 "$ENV/bin/python" -s -c "import webview, osgeo.gdal, yaml, geopandas, esda, ee" 2>/dev/null
	}

	if [ "${EHAB_REINSTALL:-0}" != "1" ] && env_ok; then
		echo "The environment is already installed in $ENV."
	else
		mkdir -p "$MAMBA_DIR"
		if [ ! -x "$MAMBA" ]; then
			echo "1/3 Downloading micromamba (the package installer of conda-forge)..."
			URL="https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-$PLATFORM"
			curl -fL --progress-bar -o "$MAMBA.part" "$URL" || fail "micromamba could not be downloaded."
			SHA=$(curl -fsL "$URL.sha256") || fail "the checksum of micromamba could not be downloaded."
			[ "$(shasum -a 256 "$MAMBA.part" | cut -d' ' -f1)" = "$SHA" ] || fail "the download of micromamba is damaged."
			chmod +x "$MAMBA.part"
			mv "$MAMBA.part" "$MAMBA"
		else
			echo "1/3 micromamba is already in $MAMBA_DIR."
		fi

		echo "2/3 List of packages: $CODE/environment.yml"
		echo "3/3 Installing Python, GDAL, GRASS GIS, R, Earth Engine and MonitorEO (a few minutes)..."
		rm -rf "$ENV"
		export MAMBA_ROOT_PREFIX="$MAMBA_DIR" PYTHONNOUSERSITE=1
		"$MAMBA" create -y -p "$ENV" -f "$CODE/environment.yml" || fail "the environment could not be installed."
		"$MAMBA" clean -y --all >/dev/null 2>&1 || true  # the cache of the packages, already installed
		env_ok || fail "the environment was installed, but something is missing (Python, R, GRASS GIS or pywebview)."
	fi

	# The app of the Desktop: an .app bundle with a script that opens eHabitat+ with the Python of the environment
	echo "Creating the app \"eHabitat+\" in $DESKTOP..."
	mkdir -p "$DESKTOP"
	rm -rf "$APP"
	mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"
	ICON="$CODE/app/ehabitat_app/ui/icon.png"
	if [ -f "$ICON" ]; then
		SET="$DIR/icon.iconset"
		rm -rf "$SET" && mkdir -p "$SET"
		for s in 16 32 128 256 512; do
			sips -z $s $s "$ICON" --out "$SET/icon_${s}x${s}.png" >/dev/null 2>&1 || true
			sips -z $((s * 2)) $((s * 2)) "$ICON" --out "$SET/icon_${s}x${s}@2x.png" >/dev/null 2>&1 || true
		done
		iconutil -c icns "$SET" -o "$APP/Contents/Resources/icon.icns" >/dev/null 2>&1 || true
		rm -rf "$SET"
	fi
	cat > "$APP/Contents/Info.plist" <<-EOF
	<?xml version="1.0" encoding="UTF-8"?>
	<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
	<plist version="1.0">
	<dict>
		<key>CFBundleName</key><string>eHabitat+</string>
		<key>CFBundleDisplayName</key><string>eHabitat+</string>
		<key>CFBundleIdentifier</key><string>es.ugr.ehabitat</string>
		<key>CFBundleExecutable</key><string>ehabitat</string>
		<key>CFBundleIconFile</key><string>icon</string>
		<key>CFBundlePackageType</key><string>APPL</string>
		<key>CFBundleVersion</key><string>1</string>
		<key>NSHighResolutionCapable</key><true/>
	</dict>
	</plist>
	EOF
	cat > "$APP/Contents/MacOS/ehabitat" <<-EOF
	#!/bin/bash
	# Opens eHabitat+ in its window (made by its installer); its messages go to ehabitat.log.
	cd "$CODE"
	export PYTHONNOUSERSITE=1 EHAB_MAMBA="$MAMBA"
	exec "$ENV/bin/python" -s "$CODE/app/ehabitat.py" >> "$DIR/ehabitat.log" 2>&1
	EOF
	chmod +x "$APP/Contents/MacOS/ehabitat"

	echo
	echo "Done. To open eHabitat+: double click \"eHabitat+\" on the Desktop. You can close the Terminal now."
	if [ "${EHAB_NO_OPEN:-0}" != "1" ]; then
		open "$APP"
	fi
}

main "$@"
