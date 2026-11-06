#!/bin/bash
# Installer of the eHabitat+ app for Ubuntu (LTS) and other Linux distributions on x86_64 (see the page Desktop app
# of the wiki of eHabpy). Run it from the folder app/installer/linux of a copy of the repository:
#
#   bash app/installer/linux/install.sh
#
# or, without downloading anything first, paste this in a terminal (wget comes with Ubuntu; curl works too):
#
#   bash -c "$(wget -qO- https://raw.githubusercontent.com/javimarlop/eHabpy/master/app/installer/linux/install.sh)"
#
# It installs in ~/eHabitat, without administrator rights (no sudo, no apt):
#   - micromamba (the package installer of conda-forge), in ~/eHabitat/micromamba;
#   - the environment of eHabitat+ (Python, GDAL, GRASS GIS, R and its packages, Earth Engine, MonitorEO and
#     pywebview, from the environment.yml of the repository, plus Qt WebEngine for the window), in ~/eHabitat/env;
#   - the code of eHabitat+ in ~/eHabitat/eHabpy, only if the installer is not run from a copy of the repository
#     (if it is, the app uses that copy as its working folder: inVars, segm_optim/pa, results...).
# It adds "eHabitat+" to the applications menu (and to the Desktop, if there is one), which opens the app in its
# window (its messages go to ~/eHabitat/ehabitat.log), and opens it. If the environment already exists and is
# complete, it only creates the launchers; to install it again, delete ~/eHabitat/env or run the installer with
# EHAB_REINSTALL=1. The app keeps itself updated from GitHub every time it opens (see Updates in that page of the wiki).
# Other variables: EHAB_HOME (installation folder), EHAB_BRANCH (branch of GitHub of the code it downloads) and
# EHAB_NO_OPEN=1 (do not open the app at the end). To uninstall: rm -rf ~/eHabitat and remove ehabitat.desktop from
# ~/.local/share/applications and from the Desktop.
# Everything is inside main, which is called at the end: a broken download runs nothing.

main() {
	set -euo pipefail
	REPO=javimarlop/eHabpy
	BRANCH="${EHAB_BRANCH:-master}"
	DIR="${EHAB_HOME:-$HOME/eHabitat}"
	ENV="$DIR/env"
	MAMBA_DIR="$DIR/micromamba"
	MAMBA="$MAMBA_DIR/micromamba"
	# the window of the app (pywebview) on Linux: Qt 5 WebEngine of conda-forge, inside the environment
	GUI_PACKAGES="pyqtwebengine qtpy"

	fail() {
		echo
		echo "ERROR: $1"
		echo "Check the internet connection and run the installer again. The messages above say what failed."
		exit 1
	}
	trap 'fail "the installation stopped at line $LINENO."' ERR

	fetch() {  # URL FILE: download with curl or, without it, wget
		if command -v curl >/dev/null 2>&1; then
			curl -fL --progress-bar -o "$2" "$1"
		elif command -v wget >/dev/null 2>&1; then
			wget -q --show-progress -O "$2" "$1"
		else
			fail "neither curl nor wget is installed (sudo apt install wget)."
		fi
	}
	fetch_text() {
		if command -v curl >/dev/null 2>&1; then curl -fsSL "$1"; else wget -qO- "$1"; fi
	}

	echo "Installer of eHabitat+"
	echo "======================"
	echo
	echo "It installs in $DIR, without administrator rights, what eHabitat+ needs: Python, GDAL, GRASS GIS, R,"
	echo "Google Earth Engine, MonitorEO and Qt for its window (about 1 GB to download and 4 GB on disk, only the"
	echo "first time)."
	echo

	case "$(uname -s)-$(uname -m)" in
		Linux-x86_64) PLATFORM=linux-64 ;;
		Linux-*) fail "GRASS GIS of conda-forge is only built for x86_64 Linux, and this computer is $(uname -m)." ;;
		*) fail "this installer is for Linux; use that of app/installer/macos or app/installer/windows." ;;
	esac
	mkdir -p "$DIR"

	# Code: the copy of the repository that contains this installer or, if none, that of ~/eHabitat/eHabpy
	HERE=""
	if [ -n "${BASH_SOURCE[0]:-}" ] && [ -f "${BASH_SOURCE[0]}" ]; then
		HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
	fi
	if [ -n "$HERE" ] && [ -f "$HERE/../../ehabitat.py" ] && [ -f "$HERE/../../../environment.yml" ]; then
		CODE="$(cd "$HERE/../../.." && pwd)"
	else
		CODE="$DIR/eHabpy"
		if [ -f "$CODE/app/ehabitat.py" ]; then
			echo "Code of eHabitat+: $CODE (already downloaded; left as it is)"
		else
			echo "Downloading the code of eHabitat+ (github.com/$REPO, branch $BRANCH)..."
			TMP="$DIR/download"
			rm -rf "$TMP" && mkdir -p "$TMP"
			fetch "https://github.com/$REPO/archive/refs/heads/$BRANCH.tar.gz" "$TMP/eHabpy.tar.gz" \
				|| fail "the code of eHabitat+ could not be downloaded."
			tar -xzf "$TMP/eHabpy.tar.gz" -C "$TMP" || fail "the downloaded code is damaged."
			SRC="$(find "$TMP" -maxdepth 1 -mindepth 1 -type d | head -n 1)"
			[ -f "$SRC/app/ehabitat.py" ] || fail "the branch $BRANCH of $REPO has no app (app/ehabitat.py)."
			mv "$SRC" "$CODE"
			rm -rf "$TMP"
		fi
	fi
	echo "Code of eHabitat+: $CODE"

	env_ok() {
		[ -x "$ENV/bin/python" ] && [ -x "$ENV/bin/Rscript" ] && [ -x "$ENV/bin/grass" ] \
			&& PYTHONNOUSERSITE=1 "$ENV/bin/python" -s -c \
				"import webview, osgeo.gdal, yaml, geopandas, esda, ee, qtpy, PyQt5.QtWebEngineWidgets" 2>/dev/null
	}

	if [ "${EHAB_REINSTALL:-0}" != "1" ] && env_ok; then
		echo "The environment is already installed in $ENV."
	else
		mkdir -p "$MAMBA_DIR"
		if [ ! -x "$MAMBA" ]; then
			echo "1/3 Downloading micromamba (the package installer of conda-forge)..."
			URL="https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-$PLATFORM"
			fetch "$URL" "$MAMBA.part" || fail "micromamba could not be downloaded."
			SHA=$(fetch_text "$URL.sha256") || fail "the checksum of micromamba could not be downloaded."
			[ "$(sha256sum "$MAMBA.part" | cut -d' ' -f1)" = "$(echo "$SHA" | tr -d '[:space:]')" ] \
				|| fail "the download of micromamba is damaged."
			chmod +x "$MAMBA.part"
			mv "$MAMBA.part" "$MAMBA"
		else
			echo "1/3 micromamba is already in $MAMBA_DIR."
		fi

		echo "2/3 List of packages: $CODE/environment.yml (and $GUI_PACKAGES for the window)"
		echo "3/3 Installing Python, GDAL, GRASS GIS, R, Earth Engine, MonitorEO and Qt (a few minutes)..."
		rm -rf "$ENV"
		export MAMBA_ROOT_PREFIX="$MAMBA_DIR" PYTHONNOUSERSITE=1
		# shellcheck disable=SC2086
		"$MAMBA" create -y -p "$ENV" -f "$CODE/environment.yml" $GUI_PACKAGES || fail "the environment could not be installed."
		"$MAMBA" clean -y --all >/dev/null 2>&1 || true  # the cache of the packages, already installed
		env_ok || fail "the environment was installed, but something is missing (Python, R, GRASS GIS, Qt or pywebview)."
	fi

	# Launchers: a script that opens the app with the Python of the environment, and its entries in the menu
	echo "Creating the launchers of eHabitat+..."
	LAUNCHER="$DIR/ehabitat.sh"
	cat > "$LAUNCHER" <<-EOF
	#!/bin/bash
	# Opens eHabitat+ in its window (made by its installer); its messages go to ehabitat.log.
	cd "$CODE"
	export PYTHONNOUSERSITE=1 EHAB_MAMBA="$MAMBA" PYWEBVIEW_GUI=qt
	export QT_QPA_PLATFORM="\${QT_QPA_PLATFORM:-xcb}"  # Qt 5 of conda-forge draws through X11 (XWayland on Wayland)
	exec "$ENV/bin/python" -s "$CODE/app/ehabitat.py" "\$@" >> "$DIR/ehabitat.log" 2>&1
	EOF
	chmod +x "$LAUNCHER"

	ENTRY="$(cat <<-EOF
	[Desktop Entry]
	Type=Application
	Name=eHabitat+
	Comment=Habitat Functional Types: inputs, HFT model, MonitorEO and results
	Exec="$LAUNCHER"
	Icon=$CODE/app/ehabitat_app/ui/icon.png
	Terminal=false
	Categories=Science;Education;Geography;
	EOF
	)"
	APPS="${XDG_DATA_HOME:-$HOME/.local/share}/applications"
	mkdir -p "$APPS"
	printf '%s\n' "$ENTRY" > "$APPS/ehabitat.desktop"
	chmod +x "$APPS/ehabitat.desktop"
	update-desktop-database "$APPS" >/dev/null 2>&1 || true
	DESKTOP="$(xdg-user-dir DESKTOP 2>/dev/null || true)"
	[ -n "$DESKTOP" ] || DESKTOP="$HOME/Desktop"
	if [ -d "$DESKTOP" ] && [ "$DESKTOP" != "$HOME" ]; then
		printf '%s\n' "$ENTRY" > "$DESKTOP/ehabitat.desktop"
		chmod +x "$DESKTOP/ehabitat.desktop"
		# GNOME (Ubuntu) only opens the launchers of the Desktop marked as trusted
		gio set "$DESKTOP/ehabitat.desktop" metadata::trusted true >/dev/null 2>&1 || true
	fi

	echo
	echo "Done. To open eHabitat+: look for \"eHabitat+\" in the applications menu (or on the Desktop), or run"
	echo "  $LAUNCHER"
	echo "With --browser the app opens in the web browser instead of its window."
	if [ "${EHAB_NO_OPEN:-0}" != "1" ] && [ -n "${DISPLAY:-}${WAYLAND_DISPLAY:-}" ]; then
		nohup "$LAUNCHER" >/dev/null 2>&1 &
	fi
}

main "$@"
