/* eHabitat+ app: the interface. It calls the functions of the app (api.py) through the local server. */
'use strict';

// Utilities ==========================================================================================================
const BASE = new URL('..', location.href).href; // http://127.0.0.1:<port>/<token>/
let ROOT = ''; // the copy of eHabpy
const $ = (sel, root = document) => root.querySelector(sel);
const $$ = (sel, root = document) => Array.from(root.querySelectorAll(sel));

async function api(name, ...args) {
	let res;
	try {
		const r = await fetch(BASE + 'api/' + name, {method: 'POST', headers: {'Content-Type': 'application/json'},
			body: JSON.stringify(args)});
		res = await r.json();
	} catch (e) {
		throw new Error('The app does not answer (was it closed?): ' + e.message);
	}
	if (!res || !res.ok) throw new Error(res ? res.error : 'no answer');
	return res.result;
}

function el(tag, attrs = {}, ...children) {
	const e = document.createElement(tag);
	for (const [k, v] of Object.entries(attrs || {})) {
		if (v === null || v === undefined || v === false) continue;
		if (k === 'class') e.className = v;
		else if (k === 'style' && typeof v === 'object') Object.assign(e.style, v);
		else if (k.startsWith('on') && typeof v === 'function') e.addEventListener(k.slice(2), v);
		else if (v === true) e.setAttribute(k, '');
		else e.setAttribute(k, v);
	}
	for (const c of children.flat()) {
		if (c === null || c === undefined || c === false) continue;
		e.append(c instanceof Node ? c : document.createTextNode(String(c)));
	}
	return e;
}

// the text of a Leaflet tooltip or popup: Leaflet takes a string as HTML, and names and labels come from files (the
// vector files of the user, ecoregions, comparison areas), so they always go as a text node
const leafletText = text => el('span', {}, text);

function toast(msg, kind = 'info', ms = 4500) {
	const t = el('div', {class: 'toast ' + kind}, msg);
	$('#toasts').append(t);
	setTimeout(() => t.remove(), kind === 'err' ? 9000 : ms);
}

function fail(e) {
	console.error(e);
	toast(e.message || String(e), 'err');
}

function fmt(v, d = 2) {
	if (v === null || v === undefined || v === '') return '—';
	const n = Number(v);
	if (!Number.isFinite(n)) return String(v);
	return n.toLocaleString('en-US', {maximumFractionDigits: d});
}

function store(key, value) {
	try {
		if (value === undefined) return localStorage.getItem('ehab.' + key);
		localStorage.setItem('ehab.' + key, value);
	} catch (e) { return null; }
}

function modal({title, body, buttons = [], wide = false}) {
	return new Promise(resolve => {
		const back = el('div', {class: 'modal-back'});
		const close = v => { back.remove(); document.removeEventListener('keydown', onKey); resolve(v); };
		const onKey = e => { if (e.key === 'Escape') close(null); };
		const foot = el('footer');
		for (const b of buttons) {
			foot.append(el('button', {class: 'btn ' + (b.cls || ''), onclick: async () => {
				const v = b.value ? await b.value() : b.result;
				if (v !== undefined) close(v);
			}}, b.label));
		}
		const box = el('div', {class: 'modal' + (wide ? ' wide' : '')},
			el('header', {}, el('h2', {}, title), el('button', {class: 'btn ghost icon right', onclick: () => close(null), title: 'Close'}, '✕')),
			el('div', {class: 'body'}, body),
			buttons.length ? foot : null);
		back.append(box);
		back.addEventListener('mousedown', e => { if (e.target === back) close(null); });
		document.addEventListener('keydown', onKey);
		document.body.append(back);
		const f = box.querySelector('input:not([type=checkbox]):not([type=radio]), select');
		if (f) setTimeout(() => f.focus(), 30);
	});
}

function confirmBox(title, text, ok = 'OK') {
	return modal({title, body: el('div', {}, text), buttons: [{label: 'Cancel', result: false}, {label: ok, cls: 'primary', result: true}]})
		.then(v => !!v);
}

function askText(title, label, value = '', placeholder = '') {
	const input = el('input', {type: 'text', value, placeholder, style: {width: '100%'}});
	const p = modal({title, body: el('label', {class: 'stack', style: {gap: '6px'}}, el('span', {class: 'muted'}, label), input),
		buttons: [{label: 'Cancel', result: null}, {label: 'OK', cls: 'primary', value: () => input.value.trim()}]});
	const btn = input.closest('.modal').querySelector('footer .btn.primary');
	input.addEventListener('keydown', e => { if (e.key === 'Enter') btn.click(); });
	return p;
}

// Files: the dialog of the window or, in a web browser, the path written by the user
async function pickFiles(kind, multiple = true, title = 'Choose files') {
	const got = await api('pick_files', kind, multiple);
	if (got !== null) return got;
	const text = await askText(title, multiple ? 'Full path of each file (separated by ;)' : 'Full path of the file');
	return text ? text.split(';').map(s => s.trim()).filter(Boolean) : [];
}

async function pickFolder() {
	const got = await api('pick_folder');
	if (got !== null) return got;
	return (await askText('Choose folder', 'Full path of the folder')) || '';
}

async function pickSave(kind, name) {
	const got = await api('pick_save', kind, name);
	if (got !== null) return got;
	return (await askText('Save as', 'Full path of the file', name)) || '';
}

// Jobs: the log of a background job, read while it runs ============================================================
class JobView {
	constructor(pre, bar, opts = {}) {
		this.pre = pre; this.bar = bar; this.opts = opts;
		this.id = null; this.since = 0; this.base = 0; this.nodes = [];
	}
	start(id, title) {
		this.id = id; this.since = 0; this.base = 0; this.nodes = []; this.step = '';
		if (this.pre) this.pre.textContent = '';
		this.title = title || '';
		if (this.opts.onStart) this.opts.onStart(id);
		this.poll();
	}
	cls(line) {
		if (line.startsWith('>>>')) return 'cmd';
		if (/\b(ERROR|Error|error:|Traceback|failed)\b/.test(line)) return 'errl';
		if (/\b(WARNING|Warning)\b/.test(line)) return 'warnl';
		if (/^(Done|Connected|Authenticated|Results in|Setup done|Everything is up to date)/.test(line.trim())) return 'okl';
		return '';
	}
	render(s) {
		if (!this.pre) return;
		let local = s.start - this.base;
		if (local < 0 || local > this.nodes.length) {
			this.pre.textContent = ''; this.nodes = []; this.base = s.start; local = 0;
		}
		const nearBottom = this.pre.scrollHeight - this.pre.scrollTop - this.pre.clientHeight < 40;
		for (const n of this.nodes.splice(local)) n.remove();
		const frag = document.createDocumentFragment();
		for (const line of s.lines) {
			const c = this.cls(line);
			const n = c ? el('span', {class: c}, line + '\n') : document.createTextNode(line + '\n');
			this.nodes.push(n); frag.append(n);
		}
		this.pre.append(frag);
		if (this.nodes.length > 6000) {
			const cut = this.nodes.length - 5000;
			for (const n of this.nodes.splice(0, cut)) n.remove();
			this.base += cut;
		}
		if (nearBottom) this.pre.scrollTop = this.pre.scrollHeight;
	}
	setBar(s) {
		if (!this.bar) return;
		this.bar.textContent = '';
		const t = s.elapsed >= 60 ? `${Math.floor(s.elapsed / 60)} min ${Math.round(s.elapsed % 60)} s` : `${Math.round(s.elapsed)} s`;
		if (s.state === 'running') {
			const last = [...s.lines].reverse().find(l => l.startsWith('>>>'));
			if (last) {
				const w = last.replace('>>> ', '').split(' ');
				const prog = /python|Rscript/i.test(w[0]) ? w[1] : w[0];
				this.step = prog && !prog.startsWith('-') ? prog.split(/[\\/]/).pop() : ''; // python -c: no script
			}
			this.bar.append(el('span', {class: 'spinner'}), el('span', {}, `${s.title}${this.step ? ' · ' + this.step : ''} · ${t}`));
		} else {
			const b = {done: ['ok', 'Finished'], failed: ['err', 'Error'], cancelled: ['warn', 'Stopped']}[s.state];
			this.bar.append(el('span', {class: 'badge ' + b[0]}, b[1]), el('span', {class: 'muted'}, `${s.title} · ${t}`));
			if (s.error && s.state === 'failed') this.bar.append(el('span', {class: 'small', style: {color: 'var(--err)'}}, s.error));
		}
	}
	async poll() {
		const id = this.id;
		let s;
		try {
			s = await api('job', id, this.since);
		} catch (e) {
			if (this.bar) this.bar.textContent = e.message;
			return;
		}
		if (id !== this.id) return;
		this.render(s);
		this.since = s.total;
		this.setBar(s);
		if (s.state === 'running') {
			setTimeout(() => this.poll(), 600);
		} else {
			this.id = null;
			if (this.opts.onEnd) this.opts.onEnd(s);
		}
	}
	running() { return this.id !== null; }
}

async function stopJob(view) {
	if (view.id !== null) await api('cancel', view.id).catch(fail);
}

// Maps ===============================================================================================================
const BASEMAPS = {
	sat: {name: 'Satellite (Esri)', url: 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
		attr: 'Imagery © Esri, Maxar, Earthstar Geographics', max: 19},
	osm: {name: 'OpenStreetMap', url: 'https://tile.openstreetmap.org/{z}/{x}/{y}.png', attr: '© OpenStreetMap contributors', max: 19},
	// global mosaic of Sentinel-2 without clouds (10 m), free for non-commercial use (CC BY-NC-SA 4.0)
	s2: {name: 'Sentinel-2 cloudless 2023 (EOX)', url: 'https://tiles.maps.eox.at/wmts/1.0.0/s2cloudless-2023_3857/default/g/{z}/{y}/{x}.jpg',
		attr: '<a href="https://s2maps.eu">Sentinel-2 cloudless - https://s2maps.eu</a> by <a href="https://eox.at">EOX IT Services GmbH</a> '
			+ '(Contains modified Copernicus Sentinel data 2023)', max: 15},
	topo: {name: 'Topographic (OpenTopoMap)', url: 'https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png',
		attr: '© OpenStreetMap, SRTM | © OpenTopoMap (CC-BY-SA)', max: 17},
	pnoa: {name: 'PNOA orthophoto (Spain)', url: 'https://www.ign.es/wmts/pnoa-ma?request=GetTile&service=WMTS&VERSION=1.0.0&Layer=OI.OrthoimageCoverage&Style=default&Format=image/jpeg&TileMatrixSet=GoogleMapsCompatible&TileMatrix={z}&TileRow={y}&TileCol={x}',
		attr: 'PNOA © Instituto Geográfico Nacional', max: 20},
	light: {name: 'Light (CARTO)', url: 'https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png',
		attr: '© OpenStreetMap © CARTO', max: 20},
	none: {name: 'No basemap'},
};

function makeMap(div, key, def, select) {
	const map = L.map(div, {preferCanvas: true, zoomSnap: 0.25, zoomDelta: 0.5, zoomControl: true, attributionControl: true}).setView([37.1, -3.4], 8);
	map.attributionControl.setPrefix('<a href="https://leafletjs.com">Leaflet</a>');
	let layer = null;
	const set = name => {
		if (!BASEMAPS[name]) name = def;
		if (layer) map.removeLayer(layer);
		layer = null;
		const b = BASEMAPS[name];
		if (b.url) layer = L.tileLayer(b.url, {attribution: b.attr, maxZoom: 21, maxNativeZoom: b.max, subdomains: 'abc'}).addTo(map);
		store('basemap.' + key, name);
		if (select) select.value = name;
	};
	if (select) {
		for (const [k, b] of Object.entries(BASEMAPS)) select.append(el('option', {value: k}, b.name));
		select.addEventListener('change', () => set(select.value));
	} else {
		const Ctl = L.Control.extend({onAdd: () => {
			const s = el('select', {class: 'small', style: {padding: '3px 6px', background: 'var(--panel)'}});
			for (const [k, b] of Object.entries(BASEMAPS)) s.append(el('option', {value: k}, b.name));
			s.addEventListener('change', () => set(s.value));
			L.DomEvent.disableClickPropagation(s);
			select = s;
			return s;
		}});
		new Ctl({position: 'topright'}).addTo(map);
	}
	set(store('basemap.' + key) || def);
	return map;
}

function bboxBounds(b) { return b ? L.latLngBounds([b[1], b[0]], [b[3], b[2]]) : null; }

// Lightbox: a plot at full size, with zoom (wheel) and pan (drag) =================================================
function lightbox(items, index) {
	let i = index, scale = 1, x = 0, y = 0;
	const img = el('img', {draggable: 'false'});
	const stage = el('div', {class: 'stage'}, img);
	const title = el('span', {class: 't'});
	const file = el('span', {class: 'small', style: {opacity: .7}});
	const box = el('div', {class: 'lightbox'},
		el('div', {class: 'bar'}, title, file, el('div', {class: 'btns right'},
			el('button', {class: 'btn small', title: 'Previous', onclick: () => go(-1)}, '←'),
			el('button', {class: 'btn small', title: 'Next', onclick: () => go(1)}, '→'),
			el('button', {class: 'btn small', title: 'Zoom in', onclick: () => zoomAt(1.4)}, '+'),
			el('button', {class: 'btn small', title: 'Zoom out', onclick: () => zoomAt(1 / 1.4)}, '−'),
			el('button', {class: 'btn small', onclick: () => fitImg()}, 'Fit'),
			el('button', {class: 'btn small', onclick: () => { scale = 1; center(); apply(); }}, '100 %'),
			el('button', {class: 'btn small', onclick: () => api('open_folder', items[i].path).catch(fail)}, 'Show file'),
			el('button', {class: 'btn small', onclick: () => close()}, '✕ Close'))),
		stage);
	const apply = () => { img.style.transform = `translate(${x}px, ${y}px) scale(${scale})`; };
	const center = () => {
		x = (stage.clientWidth - img.naturalWidth * scale) / 2;
		y = (stage.clientHeight - img.naturalHeight * scale) / 2;
	};
	const fitImg = () => {
		scale = Math.min(stage.clientWidth / img.naturalWidth, stage.clientHeight / img.naturalHeight, 1);
		center(); apply();
	};
	const zoomAt = (f, cx = stage.clientWidth / 2, cy = stage.clientHeight / 2) => {
		const ns = Math.max(0.05, Math.min(scale * f, 6));
		x = cx - (cx - x) * ns / scale; y = cy - (cy - y) * ns / scale; scale = ns; apply();
	};
	const load = () => {
		const it = items[i];
		title.textContent = it.label; file.textContent = it.file;
		img.onload = fitImg;
		img.src = BASE + it.url;
	};
	const go = d => { i = (i + d + items.length) % items.length; load(); };
	const onKey = e => {
		if (e.key === 'Escape') close();
		else if (e.key === 'ArrowRight') go(1);
		else if (e.key === 'ArrowLeft') go(-1);
	};
	const close = () => { box.remove(); document.removeEventListener('keydown', onKey); };
	stage.addEventListener('wheel', e => {
		e.preventDefault();
		const r = stage.getBoundingClientRect();
		zoomAt(e.deltaY < 0 ? 1.15 : 1 / 1.15, e.clientX - r.left, e.clientY - r.top);
	}, {passive: false});
	let drag = null;
	stage.addEventListener('mousedown', e => { drag = [e.clientX - x, e.clientY - y]; stage.style.cursor = 'grabbing'; });
	window.addEventListener('mouseup', () => { drag = null; stage.style.cursor = ''; });
	stage.addEventListener('mousemove', e => { if (drag) { x = e.clientX - drag[0]; y = e.clientY - drag[1]; apply(); } });
	stage.addEventListener('dblclick', e => { const r = stage.getBoundingClientRect(); zoomAt(2, e.clientX - r.left, e.clientY - r.top); });
	document.addEventListener('keydown', onKey);
	document.body.append(box);
	load();
}

function gallery(container, plots) {
	container.textContent = '';
	if (!plots.length) {
		container.append(el('div', {class: 'empty'}, 'No plots.'));
		return;
	}
	const items = plots.map(p => ({label: p.label, file: p.file, url: p.url,
		path: ROOT + '/' + decodeURIComponent(p.url.split('?')[0].replace(/^files\//, ''))}));
	plots.forEach((p, i) => {
		container.append(el('div', {class: 'thumb', onclick: () => lightbox(items, i), title: 'Click to enlarge'},
			el('img', {src: BASE + p.url, loading: 'lazy', alt: p.label}),
			el('div', {class: 'cap'}, p.label, el('span', {class: 'f'}, p.file))));
	});
}

// Navigation =========================================================================================================
const pages = {};
function show(name) {
	$$('.nav button').forEach(b => b.classList.toggle('active', b.dataset.page === name));
	$$('.page').forEach(p => p.classList.toggle('active', p.id === 'page-' + name));
	store('page', name);
	if (pages[name] && pages[name].onShow) pages[name].onShow();
	$('main').scrollTop = 0;
}
$$('.nav button').forEach(b => b.addEventListener('click', () => show(b.dataset.page)));

// 1. Input variables =================================================================================================
pages.inputs = (() => {
	const S = {sa: null, draft: null, sets: [], active: '', folder: null, set: null, map: null,
		saLayer: null, preview: null, eco: null};
	const addView = new JobView(null, $('#add-job'), {onEnd: s => {
		$('#add-job').classList.toggle('hidden', s.state === 'done');
		if (s.state === 'done') toast('Variables added.', 'ok');
		loadSet(S.folder);
		refreshSets();
		pages.similarity.reload(); // the variables of a set or of an area of the similarity may have changed
	}});

	function ensureMap() {
		if (!S.map) S.map = makeMap('in-map', 'inputs', 'osm');
		setTimeout(() => S.map.invalidateSize(), 50);
	}

	// Study areas ---------------------------------------------------------------------------------------------
	function renderDraft(info, areas, idError) {
		const d = S.draft;
		$('#sa-file').value = d ? d.file : '';
		const lay = $('#sa-layer');
		lay.textContent = '';
		const layers = info ? info.layers : [];
		for (const l of layers) lay.append(el('option', {value: l}, l));
		lay.value = d && d.layer ? d.layer : (layers[0] || '');
		lay.disabled = layers.length < 2;
		const col = $('#sa-col');
		col.textContent = '';
		col.append(el('option', {value: 'cat'}, 'cat (one ID per polygon)'));
		for (const f of (info ? info.fields : []).filter(f => f.integer)) col.append(el('option', {value: f.name}, f.name));
		col.value = d ? d.col_id : 'cat';
		const lab = $('#sa-label');
		lab.textContent = '';
		lab.append(el('option', {value: ''}, '(none)'));
		for (const f of (info ? info.fields : []).filter(f => !f.integer)) lab.append(el('option', {value: f.name}, f.name));
		lab.value = d && d.label_col ? d.label_col : '';
		const facts = $('#sa-facts');
		facts.textContent = '';
		if (info) {
			facts.append(
				el('div', {class: 'fact'}, el('div', {class: 'k'}, 'Study areas (IDs)'), el('div', {class: 'v'}, fmt(areas.length, 0)),
					el('div', {class: 's'}, `${fmt(info.count, 0)} polygons`)),
				el('div', {class: 'fact'}, el('div', {class: 'k'}, 'Area'), el('div', {class: 'v'}, info.area_km2 !== null ? fmt(info.area_km2) + ' km²' : '—')),
				el('div', {class: 'fact'}, el('div', {class: 'k'}, 'CRS'), el('div', {class: 'v small'}, info.crs || 'no CRS')));
		}
		const names = areas.slice(0, 6).map(a => a.label ? `${a.id} ${a.label}` : String(a.id)).join(', ');
		$('#sa-areas-note').textContent = areas.length ? `Study areas: ${names}${areas.length > 6 ? ', …' : ''}. `
			+ 'Choose which ones to process in "HFT mapping".' : '';
		$('#sa-error').classList.toggle('hidden', !idError);
		$('#sa-error').textContent = idError || '';
	}

	function dirty() {
		const d = S.draft;
		const saved = S.sa;
		const same = saved && d && saved.file === d.file && (saved.layer || null) === (d.layer || null) && saved.col_id === d.col_id
			&& (saved.label_col || null) === (d.label_col || null);
		$('#sa-apply').disabled = !d || same;
		$('#sa-state').textContent = !d ? '' : same ? 'Used by the model.' : 'Changes not applied yet.';
	}

	async function loadStudyArea() {
		try {
			S.sa = await api('study_area');
		} catch (e) { fail(e); S.sa = null; }
		if (S.sa) {
			S.draft = {file: S.sa.file, layer: S.sa.info.layer, col_id: S.sa.col_id || 'cat', label_col: S.sa.label_col || null};
			renderDraft(S.sa.info, S.sa.areas, S.sa.info.id_error);
		} else {
			S.draft = null;
			renderDraft(null, [], null);
		}
		dirty();
		drawStudyArea();
	}

	async function inspect(file, layer, col, label) {
		// label: '' guesses the column of the names, null leaves them out
		const r = await api('inspect_vector', file, layer || null, col || 'cat', label === undefined ? '' : label);
		S.draft = {file, layer: r.info.layer, col_id: col || 'cat', label_col: r.label_col || null};
		renderDraft(r.info, r.areas, r.id_error);
		dirty();
		drawGeojson(await api('vector_preview', file, r.info.layer || null, col || 'cat', r.label_col || null).catch(() => null));
	}

	async function drawStudyArea() {
		ensureMap();
		drawGeojson(await api('study_geojson').catch(() => null));
	}

	function drawGeojson(gj) {
		ensureMap();
		if (S.saLayer) S.map.removeLayer(S.saLayer);
		S.saLayer = null;
		if (!gj || !gj.features || !gj.features.length) return;
		S.saLayer = L.geoJSON(gj, {style: {color: '#ff3b30', weight: 2.5, fill: true, fillOpacity: 0.06},
			onEachFeature: (f, layer) => {
				const p = f.properties || {};
				if (p.id === undefined || p.id === null) return;
				layer.bindTooltip(leafletText(`ID ${p.id}${p.label ? ' · ' + p.label : ''}`), {sticky: true});
			}}).addTo(S.map);
		S.map.fitBounds(S.saLayer.getBounds(), {padding: [20, 20]});
	}

	$('#sa-pick').addEventListener('click', async () => {
		try {
			const f = await pickFiles('vector', false, 'File of the study areas');
			if (f && f.length) await inspect(f[0], null, 'cat', '');
		} catch (e) { fail(e); }
	});
	$('#sa-layer').addEventListener('change', () => inspect(S.draft.file, $('#sa-layer').value, 'cat', '').catch(fail));
	$('#sa-col').addEventListener('change', () => inspect(S.draft.file, S.draft.layer, $('#sa-col').value, S.draft.label_col).catch(fail));
	$('#sa-label').addEventListener('change', () => inspect(S.draft.file, S.draft.layer, S.draft.col_id, $('#sa-label').value || null).catch(fail));
	$('#sa-apply').addEventListener('click', async () => {
		const d = S.draft;
		try {
			await api('set_study_area', d.file, d.layer || null, d.col_id, d.label_col || null);
			toast('Study areas applied (copied to segm_optim/pa). Choose which ones to process in "HFT mapping".', 'ok');
			await loadStudyArea();
			pages.model.refresh();
		} catch (e) { fail(e); }
	});

	// Sets of variables ---------------------------------------------------------------------------------------
	async function refreshSets(select) {
		const r = await api('sets');
		S.sets = r.sets; S.active = r.active;
		const s = $('#set-select');
		s.textContent = '';
		for (const x of r.sets) {
			const area = x.area_only ? (x.area === 'comparison area' ? ' · comparison area only' : ' · ecoregion only')
				: x.area === 'comparison area' ? ' · comparison area' : '';
			s.append(el('option', {value: x.folder}, `${x.label} · ${x.study} + ${x.ecoregion} eco${area}${x.folder === r.active ? ' · in use' : ''}`));
		}
		if (!r.sets.length) s.append(el('option', {value: ''}, 'inVars (empty)'));
		const want = select !== undefined ? select : (S.folder !== null ? S.folder : r.active);
		s.value = r.sets.some(x => x.folder === want) ? want : (r.sets.some(x => x.folder === r.active) ? r.active : (r.sets[0] || {}).folder || '');
		renderActive();
		return s.value;
	}

	function renderActive() {
		const box = $('#set-active');
		box.textContent = '';
		const x = S.sets.find(y => y.folder === S.folder);
		const note = $('#set-area-only');
		note.textContent = '';
		note.classList.toggle('hidden', !(x && x.area_only));
		if (x && x.area_only) {
			// only an ecoregion or comparison area: listed, but not a set of variables for the HFT model
			const what = x.area === 'comparison area' ? `the comparison area ${x.area_name || ''}`.trim() : 'an ecoregion';
			const why = `This set has only the variables of ${what} (folder ecoregion), without GeoTIFFs of a study area: `
				+ 'the HFT model cannot use it. Compare the HFTs of a scenario with it in "Similarity" (Compare with).';
			box.append(el('button', {class: 'btn small', disabled: true, title: why}, 'Use this set in the model'));
			note.append(el('div', {}, el('b', {}, 'Not a set for the HFT model. '), why));
		} else if (S.folder === S.active) {
			box.append(el('span', {class: 'badge ok'}, 'Used by the model (INVARS_FOLDER)'));
		} else {
			box.append(el('button', {class: 'btn small primary', onclick: async () => {
				try {
					await api('use_set', S.folder);
					toast(`The model will use the variables of inVars${S.folder ? '/' + S.folder : ''}.`, 'ok');
					await refreshSets(S.folder);
					pages.model.refresh();
				} catch (e) { fail(e); }
			}}, 'Use this set in the model'));
		}
	}

	function short(crs) { return crs ? crs.replace(/ \(.*\)$/, '') : '—'; }

	function table(tbl, list, scale) { // scale: study, ecoregion or null (ecoregion/study_area: no left out files)
		tbl.textContent = '';
		tbl.append(el('thead', {}, el('tr', {}, el('th', {}, ''), el('th', {}, 'Variable'), el('th', {class: 'num'}, 'Pixel'),
			el('th', {}, 'Size'), el('th', {}, 'NoData'), el('th', {}, 'CRS'), el('th', {}, ''))));
		const body = el('tbody');
		if (!list.length) body.append(el('tr', {}, el('td', {colspan: 7, class: 'empty'}, 'No variables.')));
		for (const r of list) {
			const cb = scale ? el('input', {type: 'checkbox', title: r.included ? 'Leave out of the set' : 'Include again'}) : '';
			if (scale) cb.checked = r.included;
			if (scale) cb.addEventListener('change', async () => {
				try {
					await api('set_included', S.folder, scale, r.file, cb.checked);
					await loadSet(S.folder);
					refreshSets();
				} catch (e) { fail(e); cb.checked = !cb.checked; }
			});
			body.append(el('tr', {class: r.included ? '' : 'off'},
				el('td', {}, cb),
				el('td', {}, el('div', {class: 'vname'}, r.name), r.error ? el('div', {class: 'small', style: {color: 'var(--err)'}}, r.error) : null),
				el('td', {class: 'num nowrap'}, r.res ? fmt(r.res[0]) + ' m' : '—'),
				el('td', {class: 'nowrap'}, r.width ? `${fmt(r.width, 0)} × ${fmt(r.height, 0)}` : '—'),
				el('td', {class: 'nowrap'}, r.nodata === null || r.nodata === undefined ? el('span', {class: 'badge warn'}, 'none') : fmt(r.nodata)),
				el('td', {class: 'nowrap', title: r.crs || ''}, short(r.crs)),
				el('td', {}, r.error ? '' : el('button', {class: 'btn small ghost', title: 'Show on the map', onclick: () => preview(r)}, 'Show'))));
		}
		tbl.append(body);
	}

	function checks(ul, list) {
		ul.textContent = '';
		for (const c of list || []) ul.append(el('li', {class: c.level}, c.text));
	}

	function badge(b, list, checksList) {
		const n = list.filter(r => r.included).length;
		const bad = (checksList || []).some(c => c.level === 'error');
		const warn = (checksList || []).some(c => c.level === 'warning');
		b.className = 'badge ' + (n === 0 ? '' : bad ? 'err' : warn ? 'warn' : 'ok');
		const res = list.find(r => r.included && r.res);
		b.textContent = n ? `${n} variables${res ? ' · ' + fmt(res.res[0]) + ' m' : ''}` : 'empty';
	}

	async function loadSet(folder) {
		S.folder = folder;
		try {
			S.set = await api('get_set', folder);
		} catch (e) { fail(e); return; }
		const st = S.set;
		table($('#tbl-study'), st.study, 'study');
		table($('#tbl-eco'), st.ecoregion, 'ecoregion');
		checks($('#chk-study'), st.checks.study);
		checks($('#chk-eco'), st.checks.ecoregion);
		checks($('#chk-both'), st.checks.both);
		badge($('#badge-study'), st.study, st.checks.study);
		badge($('#badge-eco'), st.ecoregion, st.checks.ecoregion);
		$('#eco-show').classList.toggle('hidden', !st.ecoregions_file);
		renderArea(st);
		renderActive();
	}

	// the area of the ecoregion grid (RESOLVE ecoregions or a comparison area) and its study_area grid
	function renderArea(st) {
		const a = st.area;
		const comp = a && a.kind === 'comparison area';
		const box = $('#eco-area');
		box.textContent = '';
		if (a && a.kind) {
			box.append(el('span', {class: 'badge ' + (comp ? 'brand' : 'info')}, comp ? 'Comparison area' : 'RESOLVE ecoregions'), ' ',
				el('b', {}, a.name || ''), a.features > 1 ? el('span', {class: 'muted'}, ` (${a.features} polygons)`) : '',
				comp && a.file ? el('div', {class: 'faint mono', title: a.file, style: {wordBreak: 'break-all'}},
					a.file + (a.layer ? ` (layer ${a.layer})` : '')) : '');
		} else if (st.ecoregion.some(r => r.included)) {
			box.append(el('span', {class: 'muted'}, 'GeoTIFFs of your own: the area they cover is not recorded (no gee_vars.json).'));
		}
		$('#eco-show').textContent = comp ? 'Show the comparison area on the map' : 'Show the ecoregions on the map';
		const es = st.eco_study || [];
		$('#eco-study').classList.toggle('hidden', !es.length);
		if (es.length) {
			const res = es.find(r => r.res);
			$('#eco-study-sum').textContent = `The study area at the pixel of this grid (ecoregion/study_area): ${es.length} variables`
				+ (res ? ` · ${fmt(res.res[0])} m` : '');
			table($('#tbl-eco-study'), es, null);
			checks($('#chk-eco-study'), st.checks.eco_study);
		}
	}

	async function preview(r) {
		ensureMap();
		const leg = $('#in-legend');
		leg.textContent = 'Loading the preview…';
		try {
			const p = await api('raster_preview', r.path);
			if (S.preview) S.map.removeLayer(S.preview);
			const b = bboxBounds(p.bbox);
			S.preview = L.imageOverlay(p.image, b, {opacity: 0.85}).addTo(S.map);
			if (S.saLayer) S.saLayer.bringToFront();
			S.map.fitBounds(b, {padding: [10, 10]});
			leg.textContent = '';
			leg.append(el('b', {}, r.name), el('span', {class: 'muted'}, fmt(p.min, 3)),
				el('div', {class: 'ramp', style: {width: '160px', background: `linear-gradient(90deg, ${p.ramp.join(',')})`}}),
				el('span', {class: 'muted'}, fmt(p.max, 3)), el('span', {class: 'faint'}, '(percentiles 2–98)'),
				el('button', {class: 'btn small ghost', onclick: () => { S.map.removeLayer(S.preview); S.preview = null; leg.textContent = ''; }}, 'Remove'));
		} catch (e) { leg.textContent = ''; fail(e); }
	}

	$('#set-select').addEventListener('change', e => loadSet(e.target.value));
	$('#set-refresh').addEventListener('click', async () => { await loadSet(await refreshSets(S.folder)); });
	$('#set-open').addEventListener('click', () => S.set && api('open_folder', S.set.path).catch(fail));
	$('#set-new').addEventListener('click', async () => {
		const name = await askText('New set of variables', 'Name of the folder in inVars (letters, digits, - and _)');
		if (!name) return;
		try {
			const f = await api('create_set', name);
			await loadSet(await refreshSets(f));
			toast(`Created the folder inVars/${f} (and its subfolder ecoregion).`, 'ok');
		} catch (e) { fail(e); }
	});
	$$('[data-add]').forEach(b => b.addEventListener('click', async () => {
		if (addView.running()) return toast('Wait until the variables are added.', 'warn');
		const scale = b.dataset.add;
		try {
			const files = await pickFiles('tif', true, 'GeoTIFFs of the variables');
			if (!files.length) return;
			$('#add-job').classList.remove('hidden');
			addView.start(await api('add_rasters', S.folder, scale, files), 'Add variables');
		} catch (e) { fail(e); }
	}));
	$('#eco-show').addEventListener('click', async () => {
		ensureMap();
		try {
			const gj = await api('vector_geojson', S.set.ecoregions_file);
			if (S.eco) S.map.removeLayer(S.eco);
			const comp = S.set.area && S.set.area.kind === 'comparison area';
			S.eco = L.geoJSON(gj, {style: {color: comp ? '#c2187a' : '#ffcc00', weight: 2, dashArray: '6 4', fill: false}}).addTo(S.map);
			S.map.fitBounds(S.eco.getBounds());
		} catch (e) { fail(e); }
	});

	return {
		async init() {
			ensureMap();
			await loadStudyArea();
			await loadSet(await refreshSets());
		},
		onShow() { if (S.map) setTimeout(() => S.map.invalidateSize(), 50); },
		async refresh() { await loadSet(await refreshSets(S.folder)); },
		studyArea: () => S.sa,
	};
})();

// 2. Download from MonitorEO =========================================================================================
pages.gee = (() => {
	const S = {recipe: null, catalog: null, community: undefined, mo: null, areaInfo: null, areaLayer: null, check: null, checkSeq: 0};
	const view = new JobView($('#g-log'), $('#g-bar'), {
		onStart: () => { $('#g-stop').classList.remove('hidden'); $('#dot-gee').classList.remove('hidden'); },
		onEnd: s => {
			$('#g-stop').classList.add('hidden'); $('#dot-gee').classList.add('hidden');
			if (s.state === 'done' && view.kind === 'download') {
				toast('Download finished.', 'ok');
				pages.inputs.refresh();
				pages.model.refresh();
				scheduleCheck();
			}
			// a download (also one stopped or failed half-way) adds or changes folders of ecoregion variables: the
			// Similarity page reads its areas ("Compare with") again
			if (view.kind === 'download') pages.similarity.reload();
			// monitoreo may have changed (Update MonitorEO, the update when the app opens, or get_vars_gee.py before a
			// download or a plan): status() reads its version and, if it changed, the catalog again
			status();
			if (view.kind === 'connect' && s.state === 'done' && s.result) {
				$('#ee-eecu').textContent = s.result.eecu ? s.result.eecu.text || '' : '';
			}
		}});

	const QUALITY = [['clouds', 'Clouds and shadows'], ['water', 'Water'], ['snow', 'Snow'], ['aerosols', 'High aerosols']];
	const qbox = $('#r-quality');
	for (const [k, label] of QUALITY) {
		qbox.append(el('label', {}, label), el('select', {id: 'r-q-' + k},
			el('option', {value: 'true'}, 'Remove'), el('option', {value: 'false'}, 'Keep'),
			el('option', {value: 'auto'}, 'Auto (except in the variable that measures it)')));
	}

	function account() { return $('input[name=ee-acc]:checked').value; }
	function syncAccount() {
		const own = account() === 'own';
		$('#ee-project').disabled = !own; $('#ee-auth').disabled = !own;
		for (const id of ['#r-st-drive', '#r-eco-drive']) {
			$(id).disabled = !own || (id === '#r-eco-drive' && !$('#r-eco-on').checked);
			if (!own) $(id).checked = false;
		}
	}
	$$('input[name=ee-acc]').forEach(r => r.addEventListener('change', syncAccount));

	async function status() {
		try {
			const s = await api('gee_status');
			$('#mo-version').textContent = s.version || 'not installed';
			$('#mo-commit').textContent = s.commit ? s.commit.slice(0, 7) : s.editable ? 'development copy' : '';
			$('#mo-commit').title = s.commit ? `Commit ${s.commit} of github.com/monitoreo-obsnev/monitoreo` + (s.editable ? ` (development copy in ${s.editable})` : '') : '';
			const mo = `${s.version || ''}|${s.commit || ''}`;
			if (S.mo !== null && mo !== S.mo) {
				// another version of monitoreo was installed: its catalog, without opening the app again
				S.catalog = null;
				S.community = undefined;
				try {
					S.catalog = await api('gee_catalog');
					lcTable(S.catalog.landcover);
					ensureCommunity().then(() => { if (S.recipe) renderVars(); }, () => {});
					toast(`MonitorEO updated to ${s.version}${s.commit ? ' (' + s.commit.slice(0, 7) + ')' : ''}: the catalog has its variables.`, 'ok');
				} catch (e) { fail(e); }
			}
			S.mo = mo;
			$(`input[name=ee-acc][value=${s.account}]`).checked = true;
			if (!$('#ee-project').value) $('#ee-project').value = s.ee_project || '';
			const st = $('#ee-state');
			st.textContent = '';
			if (s.connected) st.append(el('span', {class: 'badge ok'}, s.shared ? 'Connected (shared account)' : 'Connected: ' + s.project));
			else st.append(el('span', {class: 'badge'}, 'Not connected'), ' ',
				el('span', {class: 'muted'}, 'the downloads connect by themselves'));
			syncAccount();
		} catch (e) { fail(e); }
	}

	// Recipe <-> form --------------------------------------------------------------------------------------------
	function fill(rec, source) {
		S.recipe = rec;
		$('#r-source').textContent = source ? 'from ' + source : '';
		$('#r-folder').value = rec.folder || '';
		const crs = String(rec.crs || 'auto');
		$('#r-crs-mode').value = crs.toLowerCase() === 'auto' ? 'auto' : 'code';
		$('#r-crs').value = crs.toLowerCase() === 'auto' ? '' : crs;
		$('#r-crs').classList.toggle('hidden', $('#r-crs-mode').value === 'auto');
		$('#r-st-res').value = rec.study.resolution ?? '';
		$('#r-st-buf').value = rec.study.buffer ?? '';
		$('#r-st-drive').checked = !!rec.study.drive;
		const eco = rec.ecoregion || {};
		$('#r-eco-on').checked = !!eco.resolution;
		$('#r-eco-res').value = eco.resolution || 100;
		$('#r-eco-crs').value = eco.crs || '';
		$('#r-eco-drive').checked = !!eco.drive;
		$(`input[name=r-area][value=${rec.comparison_area ? 'file' : 'resolve'}]`).checked = true;
		$('#r-area').value = rec.comparison_area || '';
		S.areaLayer = rec.comparison_area_layer || null;
		inspectArea();
		$('#r-drive-folder').value = rec.drive_folder || '';
		$('#r-drive-dir').value = rec.drive_dir || '';
		$('#r-p0').value = (rec.period || [])[0] || '';
		$('#r-p1').value = (rec.period || [])[1] || '';
		$('#r-season-on').checked = !!rec.season;
		$('#r-s0').value = rec.season ? rec.season[0] : '';
		$('#r-s1').value = rec.season ? rec.season[1] : '';
		const q = rec.quality || {};
		for (const [k] of QUALITY) $('#r-q-' + k).value = String(q[k] ?? true);
		$('#r-cs-on').checked = q.cloud_score !== null && q.cloud_score !== undefined;
		$('#r-cs').value = q.cloud_score ?? 0.6;
		$('#r-scene').value = q.scene_clouds ?? '';
		lcShow(rec);
		syncEco();
		syncAccount();
		renderVars();
		sizeNote();
	}

	function read() {
		const num = id => { const v = $(id).value.trim(); return v === '' ? null : parseFloat(v); };
		const qv = v => v === 'true' ? true : v === 'false' ? false : 'auto';
		const rec = JSON.parse(JSON.stringify(S.recipe || {}));
		rec.folder = $('#r-folder').value.trim();
		rec.crs = $('#r-crs-mode').value === 'auto' ? 'auto' : $('#r-crs').value.trim();
		rec.study = {resolution: num('#r-st-res'), buffer: num('#r-st-buf'), drive: $('#r-st-drive').checked};
		rec.ecoregion = $('#r-eco-on').checked
			? {resolution: num('#r-eco-res'), drive: $('#r-eco-drive').checked, crs: $('#r-eco-crs').value.trim() || null}
			: {resolution: null};
		const comp = $('#r-eco-on').checked && areaMode() === 'file';
		rec.comparison_area = comp ? $('#r-area').value.trim() || null : null;
		rec.comparison_area_layer = comp && $('#r-area-layer').options.length > 1 ? $('#r-area-layer').value || null : null;
		rec.drive_folder = $('#r-drive-folder').value.trim() || 'eHabitat_GEE';
		rec.drive_dir = $('#r-drive-dir').value.trim() || null;
		rec.period = [$('#r-p0').value, $('#r-p1').value];
		rec.season = $('#r-season-on').checked ? [$('#r-s0').value.trim(), $('#r-s1').value.trim()] : null;
		rec.quality = {};
		for (const [k] of QUALITY) rec.quality[k] = qv($('#r-q-' + k).value);
		rec.quality.cloud_score = $('#r-cs-on').checked ? num('#r-cs') : null;
		const sc = num('#r-scene');
		if (sc !== null) rec.quality.scene_clouds = sc;
		lcRead(rec);
		rec.variables = S.recipe.variables;
		if (!rec.folder) throw new Error('Write the destination folder in inVars.');
		if (!/^[A-Za-z0-9_-]+$/.test(rec.folder)) throw new Error('The folder can only have letters (without accents), digits, - and _.');
		if (!rec.study.resolution) throw new Error('Write the pixel size of the study area.');
		if (rec.crs !== 'auto' && !rec.crs) throw new Error('Write the code of the coordinate system (e.g. EPSG:25830).');
		if (!rec.period[0] || !rec.period[1]) throw new Error('Write the period.');
		if (rec.season && !rec.season.every(s => /^\d\d-\d\d$/.test(s))) throw new Error('Write the season as MM-DD (e.g. 06-01).');
		if (!rec.variables.length) throw new Error('Add at least one variable.');
		if (comp && !rec.comparison_area) throw new Error('Choose the vector file of the comparison area (or the RESOLVE ecoregions).');
		if (comp && !rec.ecoregion.resolution) throw new Error('Write the pixel size of the grid of the comparison area.');
		return rec;
	}

	function syncEco() {
		const on = $('#r-eco-on').checked;
		for (const id of ['#r-eco-res', '#r-eco-crs', '#r-eco-drive']) $(id).disabled = !on || (id === '#r-eco-drive' && account() !== 'own');
		syncArea();
	}
	$('#r-eco-on').addEventListener('change', syncEco);

	// Area of the ecoregion grid: the RESOLVE ecoregions or a comparison area (a vector file) -------------------
	function areaMode() { return $('input[name=r-area]:checked').value; }
	function syncArea() {
		const on = $('#r-eco-on').checked;
		const file = areaMode() === 'file';
		$$('input[name=r-area]').forEach(r => { r.disabled = !on; });
		$$('.r-area-f').forEach(e => e.classList.toggle('hidden', !file));
		const layers = S.areaInfo ? S.areaInfo.layers || [] : [];
		$$('.r-area-l').forEach(e => e.classList.toggle('hidden', !file || layers.length < 2));
		$('#r-area').disabled = !on; $('#r-area-pick').disabled = !on; $('#r-area-layer').disabled = !on;
		scheduleCheck();
	}
	async function inspectArea() {
		const f = $('#r-area').value.trim();
		const info = $('#r-area-info');
		const sel = $('#r-area-layer');
		S.areaInfo = null;
		info.textContent = '';
		sel.textContent = '';
		if (!f) { syncArea(); return; }
		try {
			const a = await api('gee_area', f, S.areaLayer || null);
			if (f !== $('#r-area').value.trim()) return;
			S.areaInfo = a;
			for (const l of a.layers || []) sel.append(el('option', {value: l}, l));
			if (a.layer) sel.value = a.layer;
			S.areaLayer = a.layer || null;
			info.append(el('b', {}, a.name), ` · ${fmt(a.count, 0)} ${a.count === 1 ? 'polygon' : 'polygons'} · `
				+ String(a.crs || 'no CRS').replace(/ \(.*\)$/, '') + (a.area_km2 ? ` · ${fmt(a.area_km2, 1)} km²` : ''));
		} catch (e) {
			if (f !== $('#r-area').value.trim()) return;
			info.append(el('span', {class: 'err-text'}, e.message));
		}
		syncArea();
	}
	let checkTimer = null;
	function scheduleCheck() { clearTimeout(checkTimer); checkTimer = setTimeout(() => checkArea().catch(() => {}), 350); }
	const areaText = a => !a ? '' : a.kind === 'comparison area' ? `the comparison area ${a.name || ''}`.trim()
		: a.kind === 'ecoregions' ? (a.names && a.names.length ? `the RESOLVE ecoregions (${a.names.join(', ')})` : a.name || 'the RESOLVE ecoregions')
			: 'GeoTIFFs of your own';
	// whether the download replaces the ecoregion grid of another area in the destination folder
	async function checkArea() {
		clearTimeout(checkTimer);
		const box = $('#r-area-warn');
		const folder = $('#r-folder').value.trim();
		const only = $('#g-only').value || null;
		const file = areaMode() === 'file' ? $('#r-area').value.trim() : null;
		const seq = ++S.checkSeq;
		if (!folder || !$('#r-eco-on').checked || only === 'study' || (areaMode() === 'file' && !file)) {
			S.check = null; box.classList.add('hidden');
			return null;
		}
		const layer = file && $('#r-area-layer').options.length > 1 ? $('#r-area-layer').value || null : null;
		const c = await api('gee_area_check', folder, file, layer, only);
		if (seq !== S.checkSeq) return S.check;
		S.check = c;
		box.textContent = '';
		box.classList.toggle('hidden', !c.replaces);
		if (c.replaces) {
			box.append(el('div', {class: 'stack', style: {gap: '6px'}},
				el('div', {}, el('b', {}, `inVars/${c.folder}/ecoregion already has the grid of ${areaText(c.existing)}. `),
					`Downloading ${areaText(c.new)} there replaces it (and the similarity runs compared with it could not be made again).`),
				c.suggestion ? el('div', {class: 'row'}, el('span', {}, 'To keep both, download it to another folder:'),
					el('button', {class: 'btn small', onclick: () => { $('#r-folder').value = c.suggestion; scheduleCheck(); }},
						`Use ${c.suggestion}`)) : null));
		}
		return c;
	}
	function syncOnly() {
		const eco = $('#g-only').value === 'ecoregion';
		$('#g-use').disabled = eco;
		$('#g-use-lab').title = eco ? 'Only the ecoregion or comparison area: not a set of variables for the HFT model' : '';
		$('#g-only-note').classList.toggle('hidden', !eco);
		scheduleCheck();
	}
	$$('input[name=r-area]').forEach(r => r.addEventListener('change', () => { syncArea(); if (areaMode() === 'file' && !S.areaInfo) inspectArea(); }));
	$('#r-area').addEventListener('change', () => { S.areaLayer = null; inspectArea(); });
	$('#r-area-layer').addEventListener('change', () => { S.areaLayer = $('#r-area-layer').value || null; inspectArea(); });
	$('#r-area-pick').addEventListener('click', async () => {
		try {
			const f = await pickFiles('vector', false, 'Vector file of the comparison area');
			if (!f || !f.length) return;
			$('#r-area').value = f[0];
			S.areaLayer = null;
			await inspectArea();
		} catch (e) { fail(e); }
	});
	$('#r-folder').addEventListener('input', scheduleCheck);
	$('#g-only').addEventListener('change', syncOnly);
	$('#r-crs-mode').addEventListener('change', () => $('#r-crs').classList.toggle('hidden', $('#r-crs-mode').value === 'auto'));
	$('#r-st-res').addEventListener('input', sizeNote);
	$('#r-st-buf').addEventListener('input', sizeNote);
	$('#r-drive-pick').addEventListener('click', async () => { const f = await pickFolder().catch(fail); if (f) $('#r-drive-dir').value = f; });

	function sizeNote() {
		const sa = pages.inputs.studyArea();
		const res = parseFloat($('#r-st-res').value);
		const note = $('#r-size-note');
		if (!sa || !sa.info || !sa.info.bbox || !res) { note.textContent = ''; return; }
		const [w, s, e, n] = sa.info.bbox;
		const lat = (s + n) / 2 * Math.PI / 180;
		const buf = parseFloat($('#r-st-buf').value) || 0;
		const wx = (e - w) * 111320 * Math.cos(lat) + 2 * buf, hy = (n - s) * 110540 + 2 * buf;
		const px = Math.ceil(wx / res) * Math.ceil(hy / res);
		note.textContent = `Grid of the study area: about ${fmt(Math.ceil(wx / res), 0)} × ${fmt(Math.ceil(hy / res), 0)} pixels, `
			+ `${fmt(px * 4 / 1048576, 1)} MB per variable. An ecoregion can be thousands of times larger: use a coarse pixel (100–1000 m).`;
	}

	// Land cover mask: the pixels of some land covers left out of every variable, after a land cover map of a year
	// (a condition of where of MonitorEO, see monitoreo.landcover): {landcover: 'CORINE', year: 2018, exclude:
	// ['water', 'built']}. Every map goes to the same broad land covers (the table of correspondences of MonitorEO), so
	// the land covers chosen stay when the map changes. The table comes with the catalog (landcover_maps).
	const LC = {maps: null, classes: [], rule: '', state: 'loading'};
	const conditions = where => !where ? [] : Array.isArray(where) ? where : [where];
	const isLc = c => !!c && typeof c === 'object' && !Array.isArray(c) && 'landcover' in c;
	function lcFind(key) {
		const low = s => String(s ?? '').toLowerCase().replace(/[\s_-]/g, '');
		return key ? (LC.maps || []).find(m => low(m.key) === low(key) || low(m.name) === low(key)) || null : null;
	}
	function lcYears(m) {
		if (Array.isArray(m.years)) return m.years.map(Number);
		const [a, b] = String(m.years).split('-').map(Number);  // '2016-2026'
		const out = [];
		for (let y = a; y <= b; y++) out.push(y);
		return out;
	}
	// the year of the map M closest to the end of the period of the recipe (without a period, that of MonitorEO)
	function lcYear(m) {
		const ys = lcYears(m);
		const end = parseInt(($('#r-p1').value || '').slice(0, 4), 10);
		if (!end) return m.default_year;
		return ys.reduce((best, y) => Math.abs(y - end) < Math.abs(best - end) ? y : best, ys[0]);
	}
	function lcYearOptions(m, year) {
		const sel = $('#r-lc-year');
		sel.textContent = '';
		if (!m) return;
		const ys = lcYears(m);
		for (const y of ys.slice().reverse()) sel.append(el('option', {value: y}, String(y)));
		sel.value = String(ys.includes(year) ? year : lcYear(m));
		$('#r-lc-years').textContent = `available: ${Array.isArray(m.years) ? ys.join(', ') : `every year from ${ys[0]} to ${ys[ys.length - 1]}`}`;
	}
	function lcChips() { return $$('#r-lc .chip.on').map(c => c.dataset.v); }
	function lcWarn(text) {
		$('#r-lc-warn').textContent = text || '';
		$('#r-lc-warn').classList.toggle('hidden', !text);
	}
	function lcSync() {
		const m = lcFind($('#r-lc-map').value);
		for (const e of $$('.r-lc-f')) e.classList.toggle('hidden', !m);
		if (!m) return;
		const cap = s => s ? s[0].toUpperCase() + s.slice(1) : '';
		$('#r-lc-about').textContent = `${cap(m.how)}; ${m.resolution_m} m, ${m.coverage}. Its classes in the land covers of the mask: ${m.crosswalk}.`;
		for (const chip of $$('#r-lc .chip')) {
			const c = LC.classes.find(k => k.key === chip.dataset.v) || {};
			chip.title = `${cap(c.help || '')}.\n${m.name}: ${(m.classes[chip.dataset.v] || []).join(', ') || 'none of its classes'}`;
		}
		const on = lcChips();
		$('#r-lc-note').textContent = on.length
			? `Leaves out the classes ${on.flatMap(k => m.classes[k] || []).join(', ')} of ${m.name} ${$('#r-lc-year').value}. ${LC.rule}`
			: 'Choose the land covers to leave out.';
	}
	// the table of the land cover maps (that of the catalog; null: a version of monitoreo without it)
	function lcTable(table, failed) {
		const had = !!LC.maps;
		const keep = {map: $('#r-lc-map').value, year: Number($('#r-lc-year').value), on: lcChips()};
		LC.maps = table ? table.maps : null;
		LC.classes = table ? table.classes : [];
		LC.rule = table ? table.rule || '' : '';
		LC.state = table ? 'ok' : failed ? 'failed' : 'old';
		const sel = $('#r-lc-map');
		sel.length = 1;  // None
		for (const m of LC.maps || []) sel.append(el('option', {value: m.key}, `${m.name} (${m.resolution_m} m, ${m.coverage})`));
		sel.disabled = !LC.maps;
		$('#r-lc').textContent = '';
		for (const c of LC.classes) {
			const chip = el('span', {class: 'chip', 'data-v': c.key}, c.name);
			chip.addEventListener('click', () => { chip.classList.toggle('on'); lcSync(); });
			$('#r-lc').append(chip);
		}
		if (!had) return lcShow(S.recipe || {});  // the mask of the recipe, now that the maps are known
		// another version of monitoreo: the same choice, if its maps still have it
		sel.value = lcFind(keep.map) ? lcFind(keep.map).key : '';
		lcYearOptions(lcFind(sel.value), keep.year);
		for (const chip of $$('#r-lc .chip')) chip.classList.toggle('on', keep.on.includes(chip.dataset.v));
		lcSync();
	}
	// the land cover mask of the recipe REC in the form: its land cover condition of where or, in a recipe of before,
	// its classes of ESA WorldCover 2021 (exclude_landcover) as the land covers of the mask that have them
	function lcShow(rec) {
		lcWarn('');
		if (!LC.maps) {
			$('#r-lc-map').value = '';
			lcSync();
			const has = conditions(rec.where).some(isLc) || (rec.exclude_landcover || []).length;
			if (LC.state === 'loading') return;
			lcWarn(LC.state === 'old'
				? `This version of monitoreo has no land cover maps: press "Update MonitorEO" to choose one.${has ? ' The land cover mask of the recipe is kept as it is.' : ''}`
				: `The land cover maps of MonitorEO could not be read (see "Update MonitorEO").${has ? ' The land cover mask of the recipe is kept as it is.' : ''}`);
			return;
		}
		let c = conditions(rec.where).find(isLc) || null;
		const old = (rec.exclude_landcover || []).map(Number);
		const wc = lcFind('WorldCover');
		if (!c && old.length && wc) {
			const ex = LC.classes.map(k => k.key).filter(k => (wc.classes[k] || []).some(t => old.includes(parseInt(t, 10))));
			c = {landcover: wc.key, year: 2021, exclude: ex};
			lcWarn(`The recipe left out the classes ${old.join(', ')} of ESA WorldCover 2021 (exclude_landcover, the mask of `
				+ `before): they are now its land covers ${ex.map(k => (LC.classes.find(x => x.key === k) || {}).name || k).join(', ')}, `
				+ `which leave out the classes below. Check them before downloading.`);
		}
		const m = c ? lcFind(c.landcover) : null;
		if (c && !m) lcWarn(`The land cover map ${c.landcover} of the recipe is not one of MonitorEO: choose another one, or none.`);
		const year = c && c.year !== undefined && c.year !== null && c.year !== '' ? Number(c.year) : null;
		if (m && year !== null && !lcYears(m).includes(year))
			lcWarn(`${m.name} has no map of ${year} (the year of the recipe): choose one of its years.`);
		$('#r-lc-map').value = m ? m.key : '';
		lcYearOptions(m, year ?? m?.default_year);
		const ex = new Set(conditions(c && c.exclude).map(x => String(x).trim().toLowerCase()));
		for (const chip of $$('#r-lc .chip')) chip.classList.toggle('on', ex.has(chip.dataset.v));
		lcSync();
	}
	// the land cover mask of the form into the recipe REC (where; without the table, the mask of the recipe as it was)
	function lcRead(rec) {
		if (!LC.maps) return;
		const m = lcFind($('#r-lc-map').value);
		const others = conditions(rec.where).filter(c => !isLc(c));  // conditions on variables: kept as they are
		let c = null;
		if (m) {
			const ex = lcChips();
			if (!ex.length) throw new Error('Choose the land covers to leave out with the land cover map (or no land cover map).');
			c = {landcover: m.key, year: Number($('#r-lc-year').value), exclude: ex};
		}
		const all = c ? [c, ...others] : others;
		rec.where = !all.length ? null : all.length === 1 && c ? c : all;
		delete rec.exclude_landcover;  // the mask of before (classes of ESA WorldCover 2021), now in where
	}
	$('#r-lc-map').addEventListener('change', () => { lcWarn(''); lcYearOptions(lcFind($('#r-lc-map').value), null); lcSync(); });
	$('#r-lc-year').addEventListener('change', () => { lcWarn(''); lcSync(); });
	lcSync();

	// Variables of the recipe -----------------------------------------------------------------------------------
	function kind(v) { return v.index ? 'catalog' : v.image ? 'image' : 'collection'; }
	function metrics(v) { const m = v.metric ?? v.metrics; return m === undefined || m === null ? [] : Array.isArray(m) ? m : [m]; }
	const TEMPORAL = {period: 'period', annual: 'annual', monthly: 'monthly', '16-day': '16 days'};

	function renderVars() {
		const tbl = $('#r-vars');
		tbl.textContent = '';
		tbl.append(el('thead', {}, el('tr', {}, ...['#', 'Variable', 'Source', 'Metrics', 'Aggregation', 'File name', ''].map(h => el('th', {}, h)))));
		const body = el('tbody');
		const vars = S.recipe.variables || [];
		if (!vars.length) body.append(el('tr', {}, el('td', {colspan: 7, class: 'empty'}, 'No variables: add them from the catalog of MonitorEO.')));
		let n = 0;
		vars.forEach((v, i) => {
			const k = kind(v);
			const ms = metrics(v);
			const files = Math.max(1, ms.length);
			const first = n + 1; n += files;
			const src = k === 'catalog' ? `${v.sensor || ''}${v.collection ? ' · ' + v.collection : ''}`
				: `${v.image || v.collection} · band ${v.band || '?'}${v.terrain ? ' · ' + v.terrain : ''}`;
			body.append(el('tr', {},
				el('td', {class: 'num faint'}, files > 1 ? `${first}–${n}` : String(first)),
				el('td', {}, el('b', {}, k === 'catalog' ? v.index : (v.name || v.band || '?')),
					el('div', {class: 'var-src'}, k === 'catalog' ? 'catalog of MonitorEO'
						: `${communityRecord(v.image || v.collection) ? 'GEE Community Catalog' : 'Earth Engine'} ${k === 'image' ? 'image' : 'collection'}`)),
				el('td', {class: 'small'}, src),
				el('td', {}, k === 'image' ? el('span', {class: 'faint'}, 'fixed image') : ms.map(m => el('span', {class: 'badge brand', style: {marginRight: '3px'}}, m))),
				el('td', {class: 'small'}, k === 'image' ? '—' : (TEMPORAL[v.temporal || 'period'] + (v.temporal && v.temporal !== 'period' ? ' · ' + (v.interannual || 'mean') : '')
					+ (Array.isArray(v.period) ? ` · ${v.period[0]} to ${v.period[1]}` : ''))),
				el('td', {class: 'small mono'}, v.name || el('span', {class: 'faint'}, 'automatic')),
				el('td', {class: 'nowrap'},
					el('button', {class: 'btn small ghost', title: 'Move up', disabled: i === 0, onclick: () => move(i, -1)}, '↑'),
					el('button', {class: 'btn small ghost', title: 'Move down', disabled: i === vars.length - 1, onclick: () => move(i, 1)}, '↓'),
					el('button', {class: 'btn small ghost', title: 'Edit', onclick: () => editVar(i)}, 'Edit'),
					el('button', {class: 'btn small ghost', title: 'Remove', onclick: () => { vars.splice(i, 1); renderVars(); }}, '✕'))));
		});
		tbl.append(body);
	}
	function move(i, d) { const v = S.recipe.variables; [v[i], v[i + d]] = [v[i + d], v[i]]; renderVars(); }

	function metricChips(selected) {
		const box = el('div', {class: 'chips'});
		const all = ['mean', 'median', 'min', 'max', 'std', 'p5', 'p10', 'p25', 'p75', 'p90', 'p95', 'sum', 'count', 'mode'];
		const extra = selected.filter(m => !all.includes(m));
		for (const m of all.concat(extra)) {
			const c = el('span', {class: 'chip' + (selected.includes(m) ? ' on' : ''), 'data-m': m}, m);
			c.addEventListener('click', () => c.classList.toggle('on'));
			box.append(c);
		}
		const custom = el('input', {type: 'text', placeholder: 'other percentile: p2.5', style: {width: '160px'}});
		custom.addEventListener('keydown', e => {
			if (e.key !== 'Enter') return;
			const m = custom.value.trim().toLowerCase();
			if (!/^p\d+(\.\d+)?$/.test(m)) return toast('Write a percentile such as p2.5 or p99.', 'warn');
			const c = el('span', {class: 'chip on', 'data-m': m}, m);
			c.addEventListener('click', () => c.classList.toggle('on'));
			box.insertBefore(c, custom); custom.value = '';
		});
		box.append(custom);
		return box;
	}
	const chosenMetrics = box => $$('.chip.on', box).map(c => c.dataset.m);

	function temporalCtl(v) {
		const t = el('select', {}, ...[['period', 'Period: the metric of all the images'], ['annual', 'Annual: the metric of each year, then between years'],
			['monthly', 'Monthly: of each month, between years and mean of the months'], ['16-day', 'Every 16 days: between years and mean']].map(([k, l]) => el('option', {value: k}, l)));
		t.value = v.temporal || 'period';
		const inter = el('select', {}, ...['mean', 'median', 'min', 'max'].map(k => el('option', {value: k}, k)));
		inter.value = v.interannual || 'mean';
		const sync = () => { inter.disabled = t.value === 'period'; };
		t.addEventListener('change', sync); sync();
		return {t, inter};
	}

	async function ensureCatalog() {
		if (!S.catalog) S.catalog = await api('gee_catalog');
		return S.catalog;
	}

	async function catalogDialog(existing) {
		let cat;
		try { cat = await ensureCatalog(); } catch (e) { return fail(e); }
		const v0 = existing || {};
		let chosen = existing ? cat.variables.find(r => r.variable === v0.index && r.sensor.toLowerCase() === String(v0.sensor || '').toLowerCase()
			&& (!v0.collection || r.collection === v0.collection)) : null;
		const search = el('input', {type: 'text', placeholder: 'Search: NDVI, temperature, water…', class: 'grow'});
		const sensor = el('select', {}, el('option', {value: ''}, 'All sensors'), ...cat.sensors.map(s => el('option', {value: s}, s)));
		const tbl = el('table', {class: 't'});
		const wrap = el('div', {class: 'table-wrap cat-table'}, tbl);
		const detail = el('div', {class: 'stack', style: {gap: '10px'}});
		const draw = () => {
			const q = search.value.trim().toLowerCase();
			tbl.textContent = '';
			tbl.append(el('thead', {}, el('tr', {}, ...['Variable', 'Name', 'Sensor', 'Collection', 'Pixel', 'Dates', 'Unit'].map(h => el('th', {}, h)))));
			const body = el('tbody');
			for (const r of cat.variables) {
				if (sensor.value && r.sensor !== sensor.value) continue;
				if (q && !`${r.variable} ${r.name} ${r.group} ${r.collection}`.toLowerCase().includes(q)) continue;
				const tr = el('tr', {class: 'click' + (chosen === r ? ' sel' : '')},
					el('td', {}, el('b', {}, r.variable)), el('td', {}, r.name, el('div', {class: 'var-src'}, r.group || '')),
					el('td', {}, r.sensor), el('td', {class: 'small mono'}, r.collection), el('td', {class: 'num'}, fmt(r.resolution_m, 0) + ' m'),
					el('td', {class: 'small nowrap'}, r.dates), el('td', {class: 'small'}, r.unit));
				tr.addEventListener('click', () => { chosen = r; draw(); drawDetail(); });
				body.append(tr);
			}
			tbl.append(body);
		};
		let mbox, tctl, name;
		const drawDetail = () => {
			detail.textContent = '';
			if (!chosen) { detail.append(el('div', {class: 'small muted'}, 'Choose a variable of the table.')); return; }
			mbox = metricChips(existing && existing.index === chosen.variable ? metrics(existing) : ['median']);
			tctl = temporalCtl(existing || {});
			name = el('input', {type: 'text', value: (existing && existing.name) || '', placeholder: `${chosen.variable}_… (automatic)`});
			detail.append(el('div', {class: 'notice info'}, el('div', {}, el('b', {}, `${chosen.variable} · ${chosen.sensor}`), ` — ${chosen.name}. `,
				el('span', {class: 'muted'}, `${chosen.collection}, ${fmt(chosen.resolution_m, 0)} m, ${chosen.dates}, ${chosen.unit}`))),
				el('div', {class: 'form-grid'},
					el('label', {}, 'Metrics (one file each)'), mbox,
					el('label', {}, 'Temporal aggregation'), tctl.t,
					el('label', {}, 'Between years'), tctl.inter,
					el('label', {}, 'File name'), name));
		};
		search.addEventListener('input', draw);
		sensor.addEventListener('change', draw);
		draw(); drawDetail();
		return modal({title: existing ? 'Edit variable of the catalog' : 'Add variables from the catalog of MonitorEO', wide: true,
			body: el('div', {class: 'stack', style: {gap: '10px'}}, el('div', {class: 'row'}, search, sensor), wrap, detail,
				el('div', {class: 'small faint'}, `Catalog of MonitorEO (${(cat.source || {}).imported || ''}), package monitoreo ${cat.monitoreo || ''}.`)),
			buttons: [{label: 'Cancel', result: null}, {label: existing ? 'Save' : 'Add', cls: 'primary', value: () => {
				if (!chosen) { toast('Choose a variable.', 'warn'); return undefined; }
				const ms = chosenMetrics(mbox);
				if (!ms.length) { toast('Choose at least one metric.', 'warn'); return undefined; }
				const several = cat.variables.filter(r => r.variable === chosen.variable && r.sensor === chosen.sensor).length > 1;
				const v = {index: chosen.variable, sensor: chosen.sensor};
				if (several) v.collection = chosen.collection;
				v.metric = ms.length === 1 ? ms[0] : ms;
				if (tctl.t.value !== 'period') { v.temporal = tctl.t.value; v.interannual = tctl.inter.value; }
				if (name.value.trim()) v.name = name.value.trim();
				for (const k of ['period', 'season', 'quality', 'exclude_landcover', 'where', 'resample', 'aggregate']) if (existing && existing[k] !== undefined) v[k] = existing[k];
				return v;
			}}]});
	}

	// The GEE Community Catalog (gee-community-catalog.org): its images and collections of images, from the copy of its
	// index in MonitorEO (no Earth Engine), read once; null with a version of monitoreo without it. A dataset goes into
	// the recipe as any data of Earth Engine, by its ID.
	const LEGACY = 'projects/earthengine-legacy/assets/';
	async function ensureCommunity() {
		if (S.community === undefined) {
			const c = await api('gee_community');
			S.community = c ? {...c, byId: new Map(c.datasets.map(d => [d.id.replace(LEGACY, ''), d]))} : null;
		}
		return S.community;
	}
	function communityRecord(asset) {
		return S.community && asset ? S.community.byId.get(String(asset).trim().replace(LEGACY, '')) || null : null;
	}
	const licenceOf = d => d.license === 'custom' ? d.license_text : d.license;
	const extLink = (url, text) => el('a', {href: url, onclick: e => { e.preventDefault(); api('open_url', url).catch(fail); }}, text || url);

	// Any image or collection of images of Earth Engine by its ID (of its data catalog, of the GEE Community Catalog or
	// any other asset): Read asks Earth Engine what it is (services.describe_asset, only metadata) and fills its bands,
	// how it is used and the scale, offset and unit of the band when its record gives them
	async function eeDialog(existing, fromCommunity = false) {
		const v0 = existing || {};
		const OWN = ['name', 'image', 'collection', 'band', 'scale', 'offset', 'unit', 'terrain', 'metric', 'metrics', 'temporal', 'interannual', 'period'];
		let info = null;  // what Earth Engine said of the ID (describe_asset), once read
		let com = null;
		try { com = await ensureCommunity(); } catch (e) { com = null; }
		const asset0 = v0.image || v0.collection || '';
		const UNIT_HINT = 'e.g. mm, degC, kg m-2 or 1 (dimensionless)';

		const source = el('select', {},
			el('option', {value: 'id'}, 'An ID of Earth Engine (its data catalog or any asset)'),
			el('option', {value: 'community', disabled: !com}, com ? `GEE Community Catalog (${com.datasets.length} datasets)` : 'GEE Community Catalog (press "Update MonitorEO")'));
		source.value = com && (fromCommunity || communityRecord(asset0)) ? 'community' : 'id';

		// the datasets of the GEE Community Catalog, by theme and words
		const theme = el('select', {}, el('option', {value: ''}, 'All the themes'),
			...(com ? com.themes.map(t => el('option', {value: t.theme}, `${t.name} (${t.datasets})`)) : []));
		const search = el('input', {type: 'search', class: 'grow', placeholder: 'Search in the title, ID, provider or tags'});
		const list = el('select', {size: 9, style: {width: '100%'}});
		const found = el('div', {class: 'small faint'});
		const drawList = () => {
			if (!com) return;
			const words = search.value.toLowerCase().split(/\s+/).filter(Boolean);
			const rows = com.datasets.filter(d => (!theme.value || d.theme === theme.value)
				&& words.every(w => [d.title, d.id, d.provider, d.tags].join(' ').toLowerCase().includes(w)));
			const keep = list.value;
			list.textContent = '';
			for (const d of rows) {
				list.append(el('option', {value: d.id, title: [d.id, d.provider, licenceOf(d) && 'Licence: ' + licenceOf(d)].filter(Boolean).join('\n')},
					`${d.title} (${d.type === 'image' ? 'image' : 'collection'})`));
			}
			list.value = keep;
			found.textContent = `${rows.length} of ${com.datasets.length} datasets (index of the catalog of ${com.date}, in MonitorEO). `
				+ 'Choose one to read its bands in Earth Engine.';
		};
		theme.addEventListener('change', drawList);
		search.addEventListener('input', drawList);
		const comBox = el('div', {class: 'stack', style: {gap: '6px'}}, el('div', {class: 'row'}, theme, search), list, found);
		const comLabel = el('label', {}, 'Dataset');

		const id = el('input', {type: 'text', class: 'grow', value: asset0, placeholder: 'e.g. COPERNICUS/DEM/GLO30 or projects/sat-io/open-datasets/…'});
		const readBtn = el('button', {class: 'btn small', title: 'Read in Earth Engine its type, bands, scale, offset and unit (only metadata: nothing is computed)'}, 'Read');
		const about = el('div', {class: 'small muted stack', style: {gap: '3px'}});
		const kindSel = el('select', {}, el('option', {value: 'image'}, 'Fixed image (or the mosaic of a collection of tiles, e.g. a DEM)'),
			el('option', {value: 'collection'}, 'Image collection with dates (aggregated in time)'));
		kindSel.value = existing ? kind(existing) : 'image';
		const bands = el('datalist', {id: 'ee-bands-' + Date.now()});
		const band = el('input', {type: 'text', value: v0.band || '', placeholder: 'e.g. DEM', list: bands.id});
		const bandNote = el('div', {class: 'small faint'});
		const scale = el('input', {type: 'number', value: v0.scale ?? '', step: 'any', placeholder: '1'});
		const offset = el('input', {type: 'number', value: v0.offset ?? '', step: 'any', placeholder: '0'});
		const unit = el('input', {type: 'text', value: v0.unit ?? '', placeholder: UNIT_HINT});
		const terrain = el('select', {}, el('option', {value: ''}, 'The value of the band'), el('option', {value: 'slope'}, 'Slope (degrees)'), el('option', {value: 'aspect'}, 'Aspect (degrees)'));
		terrain.value = v0.terrain || '';
		const mbox = metricChips(existing ? metrics(existing) : ['mean']);
		const tctl = temporalCtl(v0);
		// a collection is aggregated over the period of the recipe, or over its own (e.g. a dataset of other years)
		const ownP = el('input', {type: 'checkbox', checked: !!v0.period});
		const p0 = el('input', {type: 'date', value: (v0.period || [])[0] || $('#r-p0').value});
		const p1 = el('input', {type: 'date', value: (v0.period || [])[1] || $('#r-p1').value});
		const pHint = el('div', {class: 'small muted'});
		const name = el('input', {type: 'text', value: v0.name || '', placeholder: 'automatic'});
		const rowsImg = [el('label', {}, 'Computation'), terrain];
		const rowsCol = [el('label', {}, 'Metrics'), mbox, el('label', {}, 'Temporal aggregation'), tctl.t, el('label', {}, 'Between years'), tctl.inter,
			el('label', {}, 'Period'), el('div', {}, el('div', {class: 'row'}, el('label', {class: 'check'}, ownP, ' its own'), p0, '–', p1), pHint)];

		const periodHint = () => {
			p0.disabled = p1.disabled = !ownP.checked;
			const [a, b] = ownP.checked ? [p0.value, p1.value] : [$('#r-p0').value, $('#r-p1').value];
			const first = info && info.dates ? info.dates[0] : null, last = info && info.dates ? info.dates[1] : null;
			pHint.textContent = (ownP.checked ? '' : `That of the recipe, ${a} to ${b}. `)
				+ (first && b && first > b ? `Its first image is of ${first}, after the end of this period: give it its own period. `
					: last && a && last < a ? `Its last image is of ${last}, before the start of this period: give it its own period. ` : '');
		};
		ownP.addEventListener('change', periodHint);
		p0.addEventListener('change', periodHint);
		p1.addEventListener('change', periodHint);

		// the scale, offset and unit of the band, as the record of its data says; REPLACE: also empty those it does not give
		const applyBand = replace => {
			const b = info && (info.bands || []).find(x => x.band === band.value.trim());
			unit.placeholder = b && b.units && !b.unit ? `its record says ${b.units}: write it in UDUNITS` : UNIT_HINT;
			bandNote.textContent = b ? [b.description, b.units, b.resolution_m && `${Math.round(b.resolution_m)} m`, b.data_type,
				b.classes || b.bitmask ? 'codes, not quantities: use the metric mode (or min, max)' : null].filter(Boolean).join(' · ') : '';
			if (!b) return;
			for (const [input, val] of [[scale, b.band_scale], [offset, b.band_offset], [unit, b.unit]]) {
				if (val !== undefined && val !== null && val !== '') input.value = val;
				else if (replace) input.value = '';
			}
		};
		band.addEventListener('change', () => applyBand(true));

		const showAbout = d => {
			about.textContent = '';
			const head = [d.title && d.title !== d.asset ? d.title : null, d.source, d.provider, d.licence && `licence: ${d.licence}`].filter(Boolean).join(' · ');
			about.append(el('div', {}, el('b', {}, head || d.asset)));
			if (d.page || d.docs) about.append(el('div', {}, 'Page (units, scale and how to cite it): ', extLink(d.page || d.docs)));
			if (d.bands) {
				const n = d.bands.length;
				const dates = d.dates ? ` Dates: ${d.dates[0] || '?'} to ${d.dates[1] || 'its last image'}.` : '';
				about.append(el('div', {}, `${d.type === 'image' ? 'An image' : 'A collection of images'} with ${n} band${n === 1 ? '' : 's'}.${dates}`));
				// the notes of MonitorEO name the arguments of its functions: here, the fields below
				for (const note of d.notes || []) about.append(el('div', {}, note
					.replace(/\(unit, in UDUNITS\)/g, '(Unit, in UDUNITS)').replace(/band_scale and band_offset/g, 'Scale and offset')
					.replace(/\bunit, Scale and offset/g, 'Unit, Scale and offset')));
			}
		};
		const readAsset = async () => {
			const asset = id.value.trim();
			if (!asset) { toast('Write the ID of the image or collection.', 'warn'); return; }
			readBtn.disabled = true;
			about.textContent = `Reading ${asset} in Earth Engine… (the first time, it connects: a few seconds)`;
			try {
				const d = await api('gee_describe_asset', asset, account(), $('#ee-project').value.trim());
				if (id.value.trim() !== asset) return;  // another ID was written or chosen meanwhile
				info = d;
				id.value = d.asset;
				bands.textContent = '';
				for (const b of d.bands || []) bands.append(el('option', {value: b.band}, [b.description, b.units].filter(Boolean).join(' · ')));
				const names = (d.bands || []).map(b => b.band);
				const other = !names.includes(band.value.trim());
				if (other) band.value = (d.arguments || {}).band || names[0] || '';
				applyBand(other);
				if (!(existing && asset0 === d.asset)) kindSel.value = d.use || kindSel.value;  // a choice of before stays
				kindSel.options[1].disabled = d.type === 'image';  // an image has no dates
				if (d.type === 'image') kindSel.value = 'image';
				sync();
				showAbout(d);
			} catch (e) {
				if (id.value.trim() === asset) about.textContent = `${asset} could not be read in Earth Engine: ${e.message}`;
			} finally { readBtn.disabled = false; }
		};
		readBtn.addEventListener('click', readAsset);
		id.addEventListener('input', () => {
			info = null;
			kindSel.options[1].disabled = false;
			bands.textContent = '';
			bandNote.textContent = '';
			about.textContent = 'Read shows its bands, how it is used and what its record says.';
		});
		// a dataset of the GEE Community Catalog chosen: its ID and its record, then its bands, read in Earth Engine
		list.addEventListener('change', () => {
			const d = communityRecord(list.value);
			if (!d) return;
			id.value = d.id;
			info = null;
			for (const input of [band, scale, offset, unit]) input.value = '';
			kindSel.options[1].disabled = false;
			kindSel.value = d.type === 'image' ? 'image' : 'collection';
			sync();
			showAbout({asset: d.id, title: d.title, source: 'GEE Community Catalog', provider: d.provider, licence: licenceOf(d), docs: d.docs});
			readAsset();
		});

		const grid = el('div', {class: 'form-grid'},
			el('label', {}, 'Source'), source, comLabel, comBox,
			el('label', {}, 'Earth Engine ID'), el('div', {class: 'row', style: {flexWrap: 'nowrap'}}, id, readBtn),
			el('label', {}, ''), about,
			el('label', {}, 'Use'), kindSel,
			el('label', {}, 'Band'), el('div', {}, band, bands, bandNote),
			el('label', {}, 'Scale and offset'), el('div', {class: 'row'}, scale, offset),
			el('label', {}, 'Unit'), unit,
			...rowsImg, ...rowsCol, el('label', {}, 'File name'), name);
		const sync = () => {
			rowsImg.forEach(e => e.classList.toggle('hidden', kindSel.value !== 'image'));
			rowsCol.forEach(e => e.classList.toggle('hidden', kindSel.value !== 'collection'));
			[comLabel, comBox].forEach(e => e.classList.toggle('hidden', source.value !== 'community'));
			periodHint();
		};
		kindSel.addEventListener('change', sync);
		source.addEventListener('change', sync);
		drawList();
		if (communityRecord(asset0)) {
			list.value = communityRecord(asset0).id;
			const d = communityRecord(asset0);
			showAbout({asset: d.id, title: d.title, source: 'GEE Community Catalog', provider: d.provider, licence: licenceOf(d), docs: d.docs});
		} else about.textContent = asset0 ? 'Read shows its bands, how it is used and what its record says.' : '';
		sync();

		const title = existing ? 'Edit the variable' : fromCommunity ? 'Dataset of the GEE Community Catalog' : 'Other collection or image of Earth Engine';
		return modal({title, wide: true, body: el('div', {class: 'stack'},
			el('div', {class: 'small muted'}, 'Any image or collection of images of Earth Engine: of its data catalog '
				+ '(developers.google.com/earth-engine/datasets), e.g. WORLDCLIM/V1/BIO band bio01 with scale 0.1 (mean annual '
				+ 'temperature), or of the GEE Community Catalog (gee-community-catalog.org), shared by the community of Earth '
				+ 'Engine. Read shows its bands and what its record says; the land cover mask of the recipe applies to them too.'), grid),
			buttons: [{label: 'Cancel', result: null}, {label: existing ? 'Save' : 'Add', cls: 'primary', value: () => {
				if (!id.value.trim() || !band.value.trim()) { toast('Write the ID and the band (Read lists its bands).', 'warn'); return undefined; }
				const names = info ? (info.bands || []).map(b => b.band) : [];
				if (names.length && !names.includes(band.value.trim())) {
					toast(`${id.value.trim()} has no band ${band.value.trim()}: use one of ${names.join(', ')}.`, 'warn');
					return undefined;
				}
				const v = {};
				if (name.value.trim()) v.name = name.value.trim();
				v[kindSel.value] = id.value.trim();
				v.band = band.value.trim();
				if (scale.value !== '') v.scale = Number(scale.value);
				if (offset.value !== '') v.offset = Number(offset.value);
				if (unit.value.trim() && !(kindSel.value === 'image' && terrain.value)) v.unit = unit.value.trim();
				if (kindSel.value === 'image') { if (terrain.value) v.terrain = terrain.value; }
				else {
					const ms = chosenMetrics(mbox);
					if (!ms.length) { toast('Choose at least one metric.', 'warn'); return undefined; }
					v.metric = ms.length === 1 ? ms[0] : ms;
					if (tctl.t.value !== 'period') { v.temporal = tctl.t.value; v.interannual = tctl.inter.value; }
					if (ownP.checked) {
						if (!p0.value || !p1.value || p0.value > p1.value) { toast('Write its own period (first and last day).', 'warn'); return undefined; }
						v.period = [p0.value, p1.value];
					}
				}
				// its other options (season, quality, where, ...) stay as they were
				for (const [k, x] of Object.entries(v0)) if (!OWN.includes(k)) v[k] = x;
				return v;
			}}]});
	}

	async function editVar(i) {
		const v = S.recipe.variables[i];
		const nv = kind(v) === 'catalog' ? await catalogDialog(v) : await eeDialog(v);
		if (nv) { S.recipe.variables[i] = nv; renderVars(); }
	}

	$('#r-add-cat').addEventListener('click', async () => {
		const v = await catalogDialog(null);
		if (v) { S.recipe.variables.push(v); renderVars(); }
	});
	$('#r-add-ee').addEventListener('click', async () => {
		const v = await eeDialog(null);
		if (v) { S.recipe.variables.push(v); renderVars(); }
	});
	$('#r-add-com').addEventListener('click', async () => {
		const v = await eeDialog(null, true);
		if (v) { S.recipe.variables.push(v); renderVars(); }
	});
	$('#r-default').addEventListener('click', async () => {
		if (!await confirmBox('Default recipe', 'The form is replaced by the recipe of the repository (segm_optim/gee_vars.yaml).')) return;
		try { const r = await api('gee_default_recipe'); fill(r.recipe, r.source); } catch (e) { fail(e); }
	});
	$('#r-load').addEventListener('click', async () => {
		try {
			const f = await pickFiles('yaml', false, 'Recipe');
			if (!f.length) return;
			const r = await api('gee_load_recipe', f[0]); fill(r.recipe, r.source);
		} catch (e) { fail(e); }
	});
	$('#r-save').addEventListener('click', async () => {
		try {
			const rec = read();
			const f = await pickSave('yaml', (rec.folder || 'recipe') + '.yaml');
			if (!f) return;
			await api('gee_save_recipe', rec, f);
			toast('Recipe saved in ' + f, 'ok');
		} catch (e) { fail(e); }
	});

	// Connection, update, download -------------------------------------------------------------------------------
	const busy = () => view.running() && (toast('Wait until the running job finishes.', 'warn'), true);
	$('#ee-connect').addEventListener('click', async () => {
		if (busy()) return;
		try {
			view.kind = 'connect';
			view.start(await api('gee_connect', account(), $('#ee-project').value.trim()), 'Connection');
		} catch (e) { fail(e); }
	});
	$('#ee-auth').addEventListener('click', async () => {
		if (busy()) return;
		try { view.kind = 'auth'; view.start(await api('gee_auth'), 'Authentication'); } catch (e) { fail(e); }
	});
	$('#mo-update').addEventListener('click', async () => {
		if (busy()) return;
		try { view.kind = 'update'; view.start(await api('gee_update'), 'Update MonitorEO'); } catch (e) { fail(e); }
	});

	async function run(dry) {
		if (busy()) return;
		let rec;
		try { rec = read(); } catch (e) { return toast(e.message, 'warn'); }
		if (!pages.inputs.studyArea()) return toast('Choose the study area first in "Input variables".', 'warn');
		if (account() === 'own' && !$('#ee-project').value.trim()) return toast('Write your Earth Engine project or use the shared account.', 'warn');
		const opts = {dry_run: dry, only: $('#g-only').value || null, overwrite: $('#g-over').checked,
			use: $('#g-use').checked && $('#g-only').value !== 'ecoregion'};
		if (!dry) {
			// the ecoregion grid of another area in the destination folder: replace it or download to another folder
			let c = null;
			try { c = await checkArea(); } catch (e) { /* the download says it */ }
			if (c && c.replaces) {
				const v = await modal({title: 'Replace the ecoregion grid?', body: el('div', {class: 'stack', style: {gap: '8px'}},
					el('div', {}, `inVars/${c.folder}/ecoregion already has the grid of ${areaText(c.existing)}. Downloading `
						+ `${areaText(c.new)} there replaces its GeoTIFFs and its polygons.`),
					c.suggestion ? el('div', {}, `To keep both, download it to inVars/${c.suggestion}`
						+ (opts.only === 'ecoregion' ? ' (a set with only the area, for "Compare with" in "Similarity").' : '.')) : null),
				buttons: [{label: 'Cancel', result: null}, c.suggestion ? {label: `Download to ${c.suggestion}`, result: 'other'} : null,
					{label: 'Replace', cls: 'danger', result: 'replace'}].filter(Boolean)});
				if (!v) return;
				if (v === 'other') {
					$('#r-folder').value = c.suggestion;
					rec.folder = c.suggestion;
					scheduleCheck();
				}
			}
		}
		try {
			await api('gee_settings', account(), $('#ee-project').value.trim());
			view.kind = dry ? 'plan' : 'download';
			view.start(await api('gee_run', rec, opts), dry ? 'Plan' : 'Download');
		} catch (e) { fail(e); }
	}
	$('#g-plan').addEventListener('click', () => run(true));
	$('#g-run').addEventListener('click', () => run(false));
	$('#g-stop').addEventListener('click', () => stopJob(view));

	return {
		async init() {
			await status();
			try { const r = await api('gee_recipe'); fill(r.recipe, r.source); } catch (e) { fail(e); }
			syncOnly();
			// the land cover maps of the mask come with the catalog of MonitorEO
			ensureCatalog().then(c => lcTable(c.landcover), () => lcTable(null, true));
			ensureCommunity().then(() => { if (S.recipe) renderVars(); }, () => {});
		},
		attach(job) {
			const kinds = {'gee': 'download', 'gee-update': 'update', 'gee-connect': 'connect', 'gee-auth': 'auth'};
			view.kind = kinds[job.kind] || job.kind;
			view.start(job.id, job.title);
		},
		onShow() { sizeNote(); scheduleCheck(); },
	};
})();

// 3. HFT mapping =====================================================================================================
pages.model = (() => {
	const S = {params: [], saved: {}, status: null};
	const view = new JobView($('#m-log'), $('#m-bar'), {
		onStart: () => { $('#m-stop').classList.remove('hidden'); $('#m-run').disabled = true; $('#dot-model').classList.remove('hidden'); $('#m-after').classList.add('hidden'); $('#m-areas-card').classList.add('disabled-zone'); },
		onEnd: s => {
			$('#m-stop').classList.add('hidden'); $('#m-run').disabled = false; $('#dot-model').classList.add('hidden');
			$('#m-areas-card').classList.remove('disabled-zone');
			refresh();
			if (s.state === 'done') {
				toast('Model finished.', 'ok');
				$('#m-after').classList.remove('hidden');
				pages.results.reload(s.result);
				pages.similarity.reload();
			} else if (s.state === 'failed') toast('The model failed: see the messages.', 'err');
		}});

	function ctl(p) {
		let c;
		const changed = () => markChanged();
		if (p.type === 'bool') {
			c = el('input', {type: 'checkbox', id: 'p-' + p.name});
			c.checked = !!p.value;
			c.addEventListener('change', changed);
			return el('label', {class: 'check'}, c, el('span', {class: 'small muted'}, 'on'));
		}
		if (p.type === 'select') {
			c = el('select', {id: 'p-' + p.name}, ...p.choices.map(([v, l]) => el('option', {value: v}, l)));
			c.value = String(p.value);
			c.addEventListener('change', changed);
			return c;
		}
		if (p.type === 'threshold') {
			c = el('select', {id: 'p-' + p.name}, ...p.choices.map(([v, l]) => el('option', {value: v}, l)));
			c.value = p.value;
			const fx = el('input', {type: 'number', id: 'p-THRESHOLD_FIXED', min: 0.1, max: 0.9, step: 0.1, value: p.fixed});
			const sync = () => { fx.classList.toggle('hidden', c.value !== 'fixed'); };
			c.addEventListener('change', () => { sync(); changed(); });
			fx.addEventListener('input', changed);
			sync();
			return el('div', {class: 'inline'}, c, fx);
		}
		if (p.type === 'number' || p.type === 'int') {
			c = el('input', {type: 'number', id: 'p-' + p.name, min: p.min, max: p.max, step: p.step || 'any', value: p.value});
			c.addEventListener('input', changed);
			return el('div', {class: 'inline'}, c, p.unit ? el('span', {class: 'unit'}, p.unit) : null);
		}
		c = el('input', {type: 'text', id: 'p-' + p.name, value: p.value ?? '', style: {width: '100%'}});
		c.addEventListener('input', changed);
		if (p.name === 'SCENARIO') {
			c.placeholder = S.status ? S.status.scenario : '';
			c.addEventListener('change', checkScenario);
		}
		if (p.type === 'folder') {
			c.style.flex = '1';
			return el('div', {class: 'inline', style: {width: '100%', flexWrap: 'nowrap'}}, c, el('button', {class: 'btn small', onclick: async () => {
				const f = await pickFolder().catch(fail); if (f) { c.value = f; changed(); }
			}}, '…'));
		}
		return c;
	}

	function defText(p) {
		if (p.type === 'bool') return p.default ? 'on' : 'off';
		if (p.type === 'select' || p.type === 'threshold') return (p.choices.find(c => String(c[0]) === String(p.default)) || [, p.default])[1];
		if (p.default === '') return 'empty';
		return String(p.default) + (p.unit ? ' ' + p.unit : '');
	}

	function render(data) {
		S.params = data.params;
		const box = $('#m-params');
		box.textContent = '';
		for (const [g, title] of data.groups) {
			const card = el('div', {class: 'card'}, el('div', {class: 'card-head'}, el('h2', {}, title)));
			for (const p of data.params.filter(p => p.group === g)) {
				card.append(el('div', {class: 'field', id: 'f-' + p.name},
					el('label', {for: 'p-' + p.name, title: p.name}, p.label),
					el('div', {class: 'ctl'}, ctl(p), el('span', {class: 'def'}, 'Default: ' + defText(p))),
					el('div', {class: 'help'}, p.help, ' ', el('code', {class: 'faint'}, p.name))));
			}
			box.append(card);
		}
		const np = $('#m-nproc');
		np.max = data.nproc_max; np.value = data.nproc;
		$('#m-nproc-v').textContent = np.value;
		S.saved = values();
		markChanged();
	}
	$('#m-nproc').addEventListener('input', () => { $('#m-nproc-v').textContent = $('#m-nproc').value; });

	function values() {
		const out = {};
		for (const p of S.params) {
			const c = $('#p-' + p.name);
			if (!c) continue;
			out[p.name] = p.type === 'bool' ? c.checked : p.type === 'select' && p.cast === 'int' ? Number(c.value)
				: (p.type === 'number' || p.type === 'int') ? (c.value === '' ? null : Number(c.value)) : c.value;
		}
		const fx = $('#p-THRESHOLD_FIXED');
		if (fx) out.THRESHOLD_FIXED = Number(fx.value);
		return out;
	}

	function isDefault(p, v) { return String(v) === String(p.default); }

	function markChanged() {
		const v = values();
		let unsaved = 0;
		for (const p of S.params) {
			const f = $('#f-' + p.name);
			if (!f) continue;
			f.classList.toggle('changed', !isDefault(p, v[p.name]));
			if (JSON.stringify(v[p.name]) !== JSON.stringify(S.saved[p.name])) unsaved++;
			if (p.only_if) f.classList.toggle('disabled', !Object.entries(p.only_if).every(([k, want]) => v[k] === want));
		}
		if (v.THRESHOLD_METHOD === 'fixed' && v.THRESHOLD_FIXED !== S.saved.THRESHOLD_FIXED) unsaved++;
		const nd = S.params.filter(p => !isDefault(p, v[p.name])).length;
		$('#m-changed').textContent = (unsaved ? `${unsaved} unsaved changes · ` : '') + (nd ? `${nd} options not at their default` : 'all options at their default');
	}

	async function checkScenario() {
		try {
			const r = await api('scenario_check', $('#p-SCENARIO').value.trim());
			if (r.error) toast(`Scenario "${r.name}": ${r.error}`, 'warn');
		} catch (e) { fail(e); }
	}

	function facts(st) {
		const box = $('#m-facts');
		box.textContent = '';
		const f = (k, v, s) => el('div', {class: 'fact'}, el('div', {class: 'k'}, k), el('div', {class: 'v'}, v), s ? el('div', {class: 's'}, s) : null);
		box.append(
			f('Study areas', st.selected_study_area || '—', st.study_area_file
				? `${st.study_area_file} · ${(st.ids || []).length} of ${(st.all_ids || []).length} to process` : 'choose them in Input variables'),
			f('Variables', `inVars${st.invars_folder ? '/' + st.invars_folder : ''}`, `${st.variables.length} variables in config.py${st.resolution ? ' · ' + st.resolution + ' m pixels' : ''}`),
			f('Scenario', st.scenario, 'segm_optim/' + st.results_dir));
		const setup = $('#m-setup');
		setup.textContent = '';
		if (st.setup && st.setup.length) {
			setup.append(el('div', {class: 'notice warn'}, el('div', {}, el('b', {}, 'setup.py will run before the steps'),
				' (it imports the study area and the variables into GRASS): ' + st.setup.join('; ') + '.')));
		}
		if (st.error) setup.append(el('div', {class: 'notice err'}, st.error));
	}

	function steps(st) {
		const ul = $('#m-steps');
		ul.textContent = '';
		(st.steps || []).forEach((s, i) => {
			const ok = !s.reasons.length;
			ul.append(el('li', {class: ok ? 'ok' : 'todo'}, el('span', {class: 'n'}, String(i + 1)), el('span', {}, s.label),
				el('span', {class: 'badge ' + (ok ? 'ok' : 'warn')}, ok ? 'up to date' : 'to run'),
				ok ? null : el('div', {class: 'why'}, s.reasons.join('; '))));
		});
	}

	// Study areas (IDs) to process: palist.csv, saved as they are ticked
	function renderAreas() {
		const box = $('#m-areas');
		box.textContent = '';
		const sa = S.sa;
		if (!sa || !sa.areas || !sa.areas.length) {
			box.append(el('div', {class: 'empty'}, 'Choose the file of the study areas in "Input variables".'));
			$('#m-areas-sub').textContent = '';
			return;
		}
		const q = $('#m-areas-q').value.trim().toLowerCase();
		for (const a of sa.areas) {
			if (q && !`${a.id} ${a.label || ''}`.toLowerCase().includes(q)) continue;
			const item = el('label', {class: 'area' + (S.ids.has(a.id) ? ' on' : ''), title: a.label ? `ID ${a.id} · ${a.label}` : `ID ${a.id}`},
				el('input', {type: 'checkbox'}), el('span', {class: 'id'}, a.id), el('span', {class: 'nm'}, a.label || ''),
				el('span', {class: 'ha'}, a.ha !== null && a.ha !== undefined ? fmt(a.ha, 0) + ' ha' : ''));
			const cb = item.querySelector('input');
			cb.checked = S.ids.has(a.id);
			cb.addEventListener('change', () => {
				cb.checked ? S.ids.add(a.id) : S.ids.delete(a.id);
				item.classList.toggle('on', cb.checked);
				saveIds();
			});
			box.append(item);
		}
		areasSub();
	}

	function areasSub(text) {
		const n = S.sa && S.sa.areas ? S.sa.areas.length : 0;
		$('#m-areas-sub').textContent = `${S.ids.size} of ${n} selected (palist.csv)` + (text ? ' · ' + text : '');
	}

	let saveTimer = null;
	function saveIds() {
		clearTimeout(saveTimer);
		areasSub('saving…');
		saveTimer = setTimeout(async () => {
			if (!S.ids.size) { areasSub('choose at least one'); return; }
			try {
				await api('set_ids', [...S.ids]);
				areasSub('saved');
				refresh(false);
			} catch (e) { fail(e); loadAreas(); }
		}, 500);
	}

	async function loadAreas() {
		try { S.sa = await api('study_area'); } catch (e) { S.sa = null; }
		S.ids = new Set(S.sa ? S.sa.ids : []);
		renderAreas();
	}
	S.ids = new Set();
	$('#m-areas-q').addEventListener('input', renderAreas);
	$('#m-areas-all').addEventListener('click', () => { if (!S.sa) return; S.sa.areas.forEach(a => S.ids.add(a.id)); renderAreas(); saveIds(); });
	$('#m-areas-none').addEventListener('click', () => { S.ids.clear(); renderAreas(); saveIds(); });

	async function refresh(areas = true) {
		try {
			const st = await api('status');
			S.status = st;
			facts(st);
			steps(st);
			const sc = $('#p-SCENARIO');
			if (sc) sc.placeholder = st.scenario;
		} catch (e) { fail(e); }
		if (areas) await loadAreas();
	}

	$('#m-refresh').addEventListener('click', refresh);
	$('#m-save').addEventListener('click', async () => {
		try {
			const ch = await api('save_params', values(), Number($('#m-nproc').value));
			toast(ch.length ? `Saved in config.py: ${ch.join(', ')}.` : 'There were no changes.', 'ok');
			S.saved = values(); markChanged(); refresh();
		} catch (e) { fail(e); }
	});
	$('#m-defaults').addEventListener('click', async () => {
		if (!await confirmBox('Default values', 'The form gets the default values of all the options (not saved yet).')) return;
		for (const p of S.params) {
			const c = $('#p-' + p.name);
			if (!c) continue;
			if (p.type === 'bool') c.checked = !!p.default;
			else if (p.type === 'threshold') { c.value = p.default; c.dispatchEvent(new Event('change')); }
			else c.value = p.default;
		}
		markChanged();
	});
	$('#m-run').addEventListener('click', async () => {
		if (view.running()) return;
		const v = values();
		try {
			const sc = await api('scenario_check', v.SCENARIO || '');
			if (sc.error) return toast(`Scenario "${sc.name}": ${sc.error}`, 'warn');
			if (sc.exists && sc.changes && sc.changes.length) {
				const ok = await confirmBox('The scenario already has results', el('div', {class: 'stack', style: {gap: '8px'}},
					el('div', {}, `${sc.folder} has results made with other options:`),
					el('ul', {class: 'small'}, ...sc.changes.slice(0, 12).map(c => el('li', {}, c))),
					el('div', {}, `If you go on, they are moved to ${sc.folder}_previous (replacing what is there). To keep them, write another scenario.`)),
					'Go on');
				if (!ok) return;
			}
			if (saveTimer) { clearTimeout(saveTimer); saveTimer = null; await api('set_ids', [...S.ids]); }
			const id = await api('run_model', v, Number($('#m-nproc').value), $('#m-from').value || null);
			S.saved = values(); markChanged();
			view.start(id, 'HFT mapping');
			refresh();
		} catch (e) { fail(e); }
	});
	$('#m-stop').addEventListener('click', async () => {
		if (await confirmBox('Stop the model', 'The running step stops. The segmentation can resume later where it stopped.', 'Stop')) stopJob(view);
	});
	$('#m-see').addEventListener('click', () => show('results'));

	return {
		async init() { await refresh(); try { render(await api('params')); } catch (e) { fail(e); } },
		refresh,
		attach(job) { view.start(job.id, job.title); },
		onShow() { refresh(); },
	};
})();

// 4. Results =========================================================================================================
pages.results = (() => {
	const S = {list: [], data: null, map: null, hft: null, outline: null, layer: null, hidden: new Set(), want: null};

	function ensureMap() {
		if (!S.map) S.map = makeMap('x-map', 'results', 'sat', $('#x-base'));
		setTimeout(() => S.map.invalidateSize(), 60);
	}

	function opt(sel, values, current) {
		sel.textContent = '';
		for (const [v, l] of values) sel.append(el('option', {value: v}, l));
		if (values.some(([v]) => v === current)) sel.value = current;
	}

	async function loadList(want) {
		try { S.list = await api('results_list'); } catch (e) { return fail(e); }
		const empty = !S.list.length;
		$('#x-empty').classList.toggle('hidden', !empty);
		for (const id of ['#x-map-card', '#x-cut-card', '#x-sc-card', '#x-tables-card']) $(id).classList.toggle('hidden', empty);
		if (empty) return;
		const st = await api('status').catch(() => null);
		const sa = (want && want.study_area) || $('#x-sa').value || (st && st.study_area);
		opt($('#x-sa'), S.list.map(x => [x.study_area, x.study_area]), sa);
		await onSa((want && want.scenario) || (st && st.scenario));
	}

	async function onSa(scenario) {
		const rec = S.list.find(x => x.study_area === $('#x-sa').value) || S.list[0];
		opt($('#x-sc'), rec.scenarios.map(s => [s, s + (s.endsWith('_previous') ? ' (previous)' : '')]), scenario || $('#x-sc').value);
		await loadScenario();
	}

	async function loadScenario() {
		try { S.data = await api('results', $('#x-sa').value, $('#x-sc').value); } catch (e) { return fail(e); }
		const ids = S.data.ids;
		opt($('#x-id'), ids.map(i => [i, areaName(i)]), $('#x-id').value);
		$('#x-id').classList.toggle('hidden', ids.length < 2);
		$('#x-id-lab').classList.toggle('hidden', ids.length < 2);
		$('#x-meta').textContent = S.data.metadata || 'No processing_metadata.txt';
		await onId();
	}

	function areaName(i) {
		const name = S.data && S.data.labels ? S.data.labels[String(i)] : null;
		return name ? `${i} · ${name}` : `ID ${i}`;
	}

	async function onId() {
		const rec = S.data.per_id[$('#x-id').value];
		if (!rec) return;
		styleOutline();
		const cuts = rec.cuts;
		const sug = cuts.find(c => c.suggested);
		const lab = c => c.cut === 'results' ? 'Result' : `${c.cut}${c.hfts ? ' · ' + c.hfts + ' HFTs' : ''}${c.suggested ? ' · suggested' : ''}`;
		const cur = $('#x-cut').value;
		opt($('#x-cut'), cuts.map(c => [c.cut, lab(c)]), cuts.some(c => c.cut === cur) ? cur : (sug || cuts[0] || {}).cut);
		gallery($('#x-sc-plots'), rec.plots);
		tables(rec);
		await onCut();
	}

	async function onCut() {
		const rec = S.data.per_id[$('#x-id').value];
		const c = rec && rec.cuts.find(c => c.cut === $('#x-cut').value);
		if (!c) {
			gallery($('#x-cut-plots'), []);
			return;
		}
		$('#x-cut-sub').textContent = `${c.cut === 'results' ? '' : 'cut ' + c.cut}${c.hfts ? ' · ' + c.hfts + ' HFTs' : ''}${c.suggested ? ' · suggested cut' : ''}`
			+ (rec.threshold ? ` · segmentation threshold ${rec.threshold}` : '');
		gallery($('#x-cut-plots'), c.plots);
		await drawHft(c);
	}

	async function drawHft(c) {
		ensureMap();
		$('#x-map-title').textContent = `HFT map${c.cut === 'results' ? '' : ' · ' + c.cut}${S.data.ids.length > 1 ? ' · ' + areaName($('#x-id').value) : ''}`;
		$('#x-map-sub').textContent = 'loading…';
		let h;
		try {
			h = await api('hft_layer', $('#x-sa').value, $('#x-sc').value, $('#x-id').value, c.cut);
		} catch (e) { $('#x-map-sub').textContent = ''; return fail(e); }
		S.hft = h; S.hidden = new Set();
		S.colour = Object.fromEntries(h.legend.map(l => [l.hft, l.colour]));
		const area = Object.fromEntries(h.legend.map(l => [l.hft, l]));
		if (S.layer) S.map.removeLayer(S.layer);
		S.layer = L.geoJSON(h.geojson, {
			style: f => style(f),
			onEachFeature: (f, layer) => {
				layer.on('click', e => {
					L.popup({maxWidth: 380}).setLatLng(e.latlng).setContent(popup(f.properties, area[f.properties.hft], h)).openOn(S.map);
				});
			},
		}).addTo(S.map);
		$('#x-map-sub').textContent = `${h.legend.length} HFTs · ${h.geojson.features.length} segments · ${h.colours_from === 'tree' ? 'colours of the plots' : 'Set1 colours (run the plots again for those of the tree)'}`;
		legend(h);
		if (!S.outline || S.outlineFor !== $('#x-sa').value) await outline();
		if (S.outline) S.outline.bringToFront();
		const key = $('#x-sa').value + '/' + $('#x-sc').value + '/' + $('#x-id').value;
		if (S.fittedFor !== key) { fit(); S.fittedFor = key; }
	}

	function style(f) {
		const hid = S.hidden.has(f.properties.hft);
		const op = Number($('#x-op').value) / 100;
		return {fillColor: S.colour[f.properties.hft] || '#888', fillOpacity: hid ? 0 : op,
			color: $('#x-edges').checked && !hid ? 'rgba(20,20,20,.55)' : 'transparent', weight: 0.6, opacity: 1};
	}
	function restyle() { if (S.layer) S.layer.setStyle(style); }

	function popup(p, a, h) {
		const box = el('div');
		box.append(el('div', {class: 'row'}, el('span', {style: {display: 'inline-block', width: '16px', height: '12px', borderRadius: '3px', background: a ? a.colour : '#888'}}),
			el('b', {}, 'HFT ' + p.hft), el('span', {class: 'muted small'}, a && a.ha !== null ? `${fmt(a.ha, 1)} ha · ${fmt(a.pct, 1)} % · ${a.segments} segments` : '')));
		box.append(el('div', {class: 'small muted'}, `Segment ${p.segm_id}${p.ha !== undefined ? ' · ' + fmt(p.ha, 2) + ' ha' : ''}`));
		const prof = h.profiles[String(p.hft)];
		if (prof) {
			const grid = el('div', {class: 'prof'});
			h.variables.forEach((v, i) => {
				const x = prof[i];
				grid.append(el('span', {class: 'lab', title: v}, shortName(v)),
					el('span', {class: 'bar'}, el('i', {style: {width: `${Math.round((x ?? 0) * 100)}%`, background: a ? a.colour : '#888'}})),
					el('span', {class: 'small num'}, x === null ? '—' : fmt(x, 2)));
			});
			box.append(el('div', {class: 'small faint', style: {marginTop: '6px'}}, 'Profile of the HFT (mean of each variable, rescaled 0–1 among the HFTs)'), grid);
		}
		return box;
	}

	function shortName(v) {
		const s = v.replace(/_sinAgua.*$/, '');
		return s.length > 28 ? s.slice(0, 27) + '…' : s;
	}

	function legend(h) {
		const box = $('#x-legend');
		box.textContent = '';
		box.classList.remove('hidden');
		box.append(el('h4', {}, 'HFTs (click: hide/show)'));
		for (const l of h.legend) {
			const row = el('div', {class: 'legend-row' + (S.hidden.has(l.hft) ? ' dim' : '')},
				el('span', {class: 'sw', style: {background: l.colour}}), el('b', {}, l.hft),
				el('span', {class: 'ha'}, l.ha !== null ? `${fmt(l.ha, 0)} ha · ${fmt(l.pct, 1)} %` : `${l.segments} segm.`));
			row.addEventListener('click', () => {
				S.hidden.has(l.hft) ? S.hidden.delete(l.hft) : S.hidden.add(l.hft);
				row.classList.toggle('dim'); restyle();
			});
			box.append(row);
		}
		box.append(el('div', {class: 'row', style: {marginTop: '6px'}},
			el('button', {class: 'btn small ghost', onclick: () => { S.hidden.clear(); legend(h); restyle(); }}, 'All'),
			el('button', {class: 'btn small ghost', onclick: () => { h.legend.forEach(l => S.hidden.add(l.hft)); legend(h); restyle(); }}, 'None')));
	}

	async function outline() {
		if (S.outline) { S.map.removeLayer(S.outline); S.outline = null; }
		S.outlineFor = $('#x-sa').value;
		const gj = await api('study_outline', $('#x-sa').value).catch(() => null);
		if (!gj) return;
		// the study area shown in black, the others of the file in grey (click one with results to show it)
		S.outline = L.geoJSON(gj, {style: outlineStyle, onEachFeature: (f, layer) => {
			const id = f.properties.id;
			if (id === undefined || id === null) return;
			layer.bindTooltip(leafletText(`ID ${id}${f.properties.label ? ' · ' + f.properties.label : ''}`), {sticky: true});
			layer.on('click', e => {
				if (String(id) === $('#x-id').value || !S.data.ids.includes(String(id))) return;
				L.DomEvent.stopPropagation(e);
				$('#x-id').value = String(id);
				onId();
			});
		}});
		if ($('#x-outline').checked) S.outline.addTo(S.map);
	}

	function outlineStyle(f) {
		const id = f.properties.id;
		const multi = S.data && S.data.ids.length > 1;
		if (id === undefined || id === null || !multi || String(id) === $('#x-id').value) return {color: '#111', weight: 2.4, fill: false};
		const done = S.data.ids.includes(String(id));
		return {color: done ? '#f5f5f5' : '#9aa5a0', weight: done ? 1.8 : 1.2, dashArray: '5 4', fill: false};
	}
	function styleOutline() { if (S.outline) S.outline.setStyle(outlineStyle); }

	function fit() {
		S.map.invalidateSize();
		if (S.hft) S.map.fitBounds(S.hft.bounds, {padding: [16, 16]});
	}

	function tables(rec) {
		const box = $('#x-tables');
		box.textContent = '';
		for (const t of rec.tables) {
			if (t.file.startsWith('r2_') && t.file.endsWith('_hfts.csv')) continue; // long: shown in the R² plot
			// hft_means_*.csv: without the columns that identify the scenario (the same in every row), kept in the file
			// to join the tables of several study areas
			const same = t.file.startsWith('hft_means_') ? new Set(['label', 'study_area', 'scenario', 'id', 'resolution', 'threshold']) : new Set();
			const cols = t.header.map((h, j) => j).filter(j => !same.has(t.header[j]));
			const head = cols.map(j => t.header[j]);
			const tbl = el('table', {class: 't'});
			tbl.append(el('thead', {}, el('tr', {}, ...head.map(h => el('th', {title: h}, shortName(h))))));
			const sugCol = t.header.indexOf('suggested');
			const selCol = t.header.indexOf('selected');
			const body = el('tbody');
			for (const r of t.rows) {
				const hi = (sugCol >= 0 && r[sugCol] === '1') || (selCol >= 0 && r[selCol] === '1');
				body.append(el('tr', {class: hi ? 'sel' : ''}, ...cols.map(j => r[j] ?? '').map(v => el('td', {class: /^-?[\d.]+(e-?\d+)?$/.test(v) ? 'num' : ''}, /^-?\d+\.\d{4,}/.test(v) ? fmt(v, 3) : v))));
			}
			tbl.append(body);
			box.append(el('div', {class: 'stack', style: {gap: '6px'}}, el('div', {class: 'row'}, el('h3', {}, t.label), el('span', {class: 'small faint'}, t.file)),
				el('div', {class: 'table-wrap'}, tbl)));
		}
		if (!rec.tables.length) box.append(el('div', {class: 'empty'}, 'No tables.'));
	}

	$('#x-sa').addEventListener('change', () => onSa());
	$('#x-sc').addEventListener('change', () => loadScenario());
	$('#x-id').addEventListener('change', () => onId());
	$('#x-cut').addEventListener('change', onCut);
	$('#x-op').addEventListener('input', restyle);
	$('#x-edges').addEventListener('change', restyle);
	$('#x-outline').addEventListener('change', () => { if (S.outline) $('#x-outline').checked ? S.outline.addTo(S.map) : S.map.removeLayer(S.outline); });
	$('#x-fit').addEventListener('click', fit);
	$('#x-open').addEventListener('click', () => S.data && api('open_folder', S.data.path).catch(fail));
	$('#x-refresh').addEventListener('click', () => loadList());

	return {
		async onShow() {
			ensureMap();
			if (!S.loaded || S.want) { S.loaded = true; const w = S.want; S.want = null; await loadList(w); }
		},
		reload(want) { S.want = want || {}; S.fittedFor = null; if ($('#page-results').classList.contains('active')) this.onShow(); },
	};
})();

// 5. Similarity in the ecoregion =====================================================================================
pages.similarity = (() => {
	const S = {list: [], src: null, plan: null, planError: null, planSeq: 0, params: [], defaults: {}, cfg: {}, saved: null,
		runs: [], res: null, map: null, overlay: null, eco: null, area: null, ecoFor: null, saLayer: null, saFor: null, previews: new Map(),
		mapSeq: 0, fitted: null, defaultEco: null, defaultEcoPath: null, loaded: false, stale: false, full: null, fullKey: null,
		detail: null, detailSeq: 0, detailTimer: null, resFollow: null, resNative: null, resKey: ''};
	const TITLE = 'Similarity in the ecoregion';
	const view = new JobView($('#s-log'), $('#s-bar'), {
		onStart: () => { $('#s-stop').classList.remove('hidden'); $('#dot-ecosim').classList.remove('hidden'); canRun(); },
		onEnd: async s => {
			$('#s-stop').classList.add('hidden'); $('#dot-ecosim').classList.add('hidden');
			canRun();
			if (s.state === 'done') {
				const r = s.result || {};
				if (r.folder) toast('Similarity computed.', 'ok');
				else toast('similarity.py finished, but no finished run (similarity_run.json) was found: see the messages.', 'warn');
				if (r.study_area && r.scenario && key(r.study_area, r.scenario) !== $('#s-sc').value
					&& [...$('#s-sc').options].some(o => o.value === key(r.study_area, r.scenario))) {
					$('#s-sc').value = key(r.study_area, r.scenario);
					await loadSources(r.folder);
				} else {
					await loadRuns(r.folder);
					refreshPlan();
				}
			} else {
				// failed or stopped: similarity.py had already removed the folder of the run it replaces
				if (s.state === 'failed') toast('The similarity failed: see the messages.', 'err');
				await loadRuns();
				refreshPlan();
			}
		}});
	const STATUS = {
		ok: ['ok', 'ok', 'computed'],
		few_pixels: ['warn', 'few pixels', 'computed with fewer pixels than 3 × the variables: less reliable'],
		too_few_pixels: ['err', 'too few pixels', 'not computed: fewer pixels than the variables + 1 (and at least 5)'],
		no_reference_similarity: ['err', 'no reference', 'no pixel of the HFT with a similarity above 0'],
		no_similar: ['warn', 'no similar pixels', 'no pixel of the rest of the ecoregion as similar as the HFT'],
	};
	const THDI = [['TerrEco', 'Terrestrial ecoregions', 'touched by the HFTs'], ['MarEco', 'Marine ecoregions', 'touched by the HFTs'],
		['NrHFT', 'HFTs', ''], ['SIH', 'SIH', 'Shannon index of the areas of the HFTs'],
		['THD', 'THD', 'Total Habitat Diversity: ecoregions × SIH'],
		['THR', 'THR', 'Total Habitat Replaceability: median HSR / terrestrial ecoregions'], ['THDI', 'THDI', 'THD / THR']];
	const HRI = [ // column, header, decimals, title
		['HFTArea_ha', 'Area (ha)', 1, 'Area of the HFT in its map'],
		['NrPixHFT', 'HFT pixels', 0, 'Pixels of the HFT on the grid of the similarity'],
		['T', 'T', 3, 'Similarity threshold: the statistic of the similarity of the HFT pixels'],
		['NrSimPixHFT', 'HFT pixels ≥ T', 0, 'Pixels of the HFT with a similarity of at least T'],
		['NrSimPixEco', 'Similar pixels', 0, 'Pixels of the rest of the ecoregion with a similarity of at least T'],
		['SimAreaEco_km2', 'Similar area (km²)', 1, 'Area of the similar pixels'],
		['NrSimLcpPatEco', 'Patches', 0, 'Patches of similar pixels'],
		['NrSimLcpPatEqLargArea', 'Patches ≥ HFT', 0, 'Patches at least as large as the HFT'],
		['LargPatArea_km2', 'Largest patch (km²)', 2, 'Area of the largest patch'],
		['HSR', 'HSR', 2, 'Habitat Similarity Ratio: pixels of the patches at least as large as the HFT (or of the largest patch) / '
			+ 'pixels of the HFT'],
		['Aggregation', 'Aggregation', 1, 'Similar pixels per patch']];

	const key = (sa, sc) => sa + '/' + sc;
	function parts() {
		const v = $('#s-sc').value || '';
		const i = v.indexOf('/');
		return i < 0 ? [null, null] : [v.slice(0, i), v.slice(i + 1)];
	}
	const base = p => p ? String(p).split(/[\\/]/).pop() : '—';
	const mem = mb => !mb ? '' : mb >= 1024 ? `~${fmt(mb / 1024, 1)} GB` : `~${fmt(mb, 0)} MB`;
	function areaName(i) {
		const name = S.src && S.src.labels ? S.src.labels[String(i)] : null;
		return name ? `${i} · ${name}` : `ID ${i}`;
	}
	function statusBadge(st) {
		const b = STATUS[st] || ['', st || '—', ''];
		return el('span', {class: 'badge ' + b[0], title: b[2]}, b[1]);
	}
	const swatch = c => el('span', {class: 'swatch', style: {background: c || '#888'}});

	// Options -----------------------------------------------------------------------------------------------------
	function defText(p) {
		if (p.type === 'select') return (p.choices.find(c => String(c[0]) === String(p.default)) || [, p.default])[1];
		if (p.name === 'ECOSIM_RESOLUTION') return 'the pixel of the ecoregion variables';
		if (p.name === 'ECOSIM_CUT') return 'the suggested one';
		if (p.name === 'ECOSIM_AREA') return 'the ecoregion of the set of variables';
		return String(p.default) + (p.unit ? ' ' + p.unit : '');
	}

	function fillFields(params) {
		S.params = params;
		S.defaults = Object.fromEntries(params.map(p => [p.name, p.default]));
		for (const p of params) {
			const f = $('#sf-' + p.name);
			if (!f) continue;
			const lab = f.querySelector('label');
			lab.textContent = p.label;
			lab.title = p.name;
			const help = f.querySelector('.help');
			help.textContent = '';
			help.append(p.help, ' ', el('code', {class: 'faint'}, p.name));
			const def = f.querySelector('.def');
			if (def) def.textContent = 'Default: ' + defText(p);
		}
		for (const [id, name] of [['#s-stat', 'ECOSIM_STAT'], ['#s-nb', 'ECOSIM_NEIGHBOURS']]) {
			const p = params.find(q => q.name === name);
			const sel = $(id);
			const cur = sel.value;
			sel.textContent = '';
			for (const [v, l] of p.choices) sel.append(el('option', {value: v}, l));
			if (cur) sel.value = cur;
		}
	}

	function setValues(v) {
		$('#s-stat').value = String(v.ECOSIM_STAT);
		$('#s-nb').value = String(v.ECOSIM_NEIGHBOURS);
		$('#s-cover').value = v.ECOSIM_MIN_COVER;
		$('#s-eco').value = v.ECOSIM_ECOREGIONS || '';
		$('#s-marine').value = v.ECOSIM_MARINE || '';
		fillRes(Number(v.ECOSIM_RESOLUTION) || 0);
		if (S.src) fillCuts(v.ECOSIM_CUT || '');
		fillAreas(v.ECOSIM_AREA || '');
		ecoPlaceholder();
	}

	// the ecoregion or comparison area: the default one (the ecoregion of the set of variables, value '') first, then
	// the other folders of ecoregion variables of inVars (similarity.comparison_areas)
	const shortFolder = f => { const p = String(f || '').split('/'); return p.length > 3 && f.startsWith('/') ? '…/' + p.slice(-3).join('/') : f; };
	function areaLabel(a) {
		const what = a.kind === 'comparison area' ? 'comparison area' : 'ecoregions';
		const name = a.name ? `${what}: ${a.name}` : what;
		const parts = a.default ? [`The ecoregion of the set (${a.folder})`, name] : [shortFolder(a.folder), name];
		if (a.res) parts.push(`${fmt(a.res)} m`);
		if (a.reference) parts.push('with study_area');
		if ((a.missing || []).length) parts.push('missing variables: ' + a.missing.join(', '));
		else if (a.error) parts.push(a.error);
		return parts.join(' · ');
	}
	function areaValue(want) { // the value of the select for a folder of ECOSIM_AREA ('' = the default area)
		const def = ((S.src && S.src.areas) || []).find(a => a.default);
		const w = String(want || '').replace(/\\/g, '/').replace(/\/+$/, '');
		return def && w === def.folder ? '' : w;
	}
	function fillAreas(want) {
		const sel = $('#s-area');
		sel.textContent = '';
		const list = (S.src && S.src.areas) || [];
		const def = list.find(a => a.default);
		sel.append(el('option', {value: ''}, def ? areaLabel(def) : 'The ecoregion of the set of variables'));
		for (const a of list) {
			if (!a.default) sel.append(el('option', {value: a.folder, title: a.folder}, areaLabel(a)));
		}
		const w = areaValue(want);
		if (w && ![...sel.options].some(o => o.value === w)) sel.append(el('option', {value: w}, `${shortFolder(w)} (not found)`));
		sel.value = w;
	}
	function chosenArea() {
		const v = $('#s-area').value;
		return ((S.src && S.src.areas) || []).find(a => (v === '' && a.default) || (v !== '' && a.folder === v)) || null;
	}

	// pixel size: the choices of the plan (1, 2, 5 and 10 times the pixel of the ecoregion variables) or another. The
	// select is always that of the plan of the area chosen (S.resNative, S.resKey): when the area changes, the area is
	// planned at its own pixel (S.resFollow: the pixel size wanted, kept if it is a multiple of the new pixel)
	const planNative = p => p && p.native && p.native.res ? Number(p.native.res) : null;
	const resKey = p => planNative(p) ? JSON.stringify([planNative(p), p.resolutions || []]) : '';
	function native() { return S.resNative; }
	function keepRes(want, p) { // WANT (ECOSIM_RESOLUTION) if it is a whole multiple of the pixel of the plan P, else 0
		const nat = planNative(p);
		const f = want && nat ? Number(want) / nat : 0;
		return f > 1 + 1e-9 && Math.abs(f - Math.round(f)) < 1e-6 ? Number(want) : 0;
	}
	function resLabel(r) {
		return `${fmt(r.res)} m${r.factor === 1 ? ' (ecoregion variables)' : ''} · ${fmt(r.width, 0)} × ${fmt(r.height, 0)} px`
			+ (r.memory_mb ? ' · ' + mem(r.memory_mb) : '');
	}
	function fillRes(want) { // want: ECOSIM_RESOLUTION (0 = the pixel of the ecoregion variables), with the plan S.plan
		const sel = $('#s-res');
		sel.textContent = '';
		const p = S.plan;
		S.resNative = planNative(p);
		S.resKey = resKey(p);
		const list = S.resNative && p.resolutions ? p.resolutions : [];
		for (const r of list) sel.append(el('option', {value: String(r.res)}, resLabel(r)));
		if (!list.length) sel.append(el('option', {value: '0'}, 'The pixel of the ecoregion variables'));
		sel.append(el('option', {value: 'other'}, 'Other…'));
		const nat = native();
		const w = want ? Number(want) : (nat || 0);
		const same = [...sel.options].find(o => o.value !== 'other' && Math.abs(Number(o.value) - w) < 1e-9);
		if (same) sel.value = same.value;
		else if (!want) sel.value = sel.options[0].value;
		else { sel.value = 'other'; $('#s-res-other').value = want; }
		$('#s-res-other').step = nat || 1;
		$('#s-res-other').min = nat || 1;
		syncRes();
	}
	function syncRes() {
		$('#s-res-other-box').classList.toggle('hidden', $('#s-res').value !== 'other');
		resCheck();
	}
	function resMetres() { // the pixel size in metres for the plan (0 = unknown: that of config.py)
		const v = $('#s-res').value;
		if (v === 'other') return $('#s-res-other').value === '' ? null : Number($('#s-res-other').value);
		return Number(v) || native() || 0;
	}
	function resConfig() { // the value of ECOSIM_RESOLUTION
		const m = resMetres();
		const nat = native();
		return m === null ? null : nat && Math.abs(m - nat) < 1e-9 ? 0 : m;
	}
	function resError() {
		if ($('#s-res').value !== 'other') return null;
		const m = resMetres();
		const nat = native();
		if (m === null || !Number.isFinite(m) || m <= 0) return 'write the pixel size in metres';
		if (nat) {
			const f = m / nat;
			if (f < 1 - 1e-9 || Math.abs(f - Math.round(f)) > 1e-6) {
				return `the pixel size must be a whole multiple of ${fmt(nat)} m, the pixel of the ecoregion variables `
					+ `(e.g. ${fmt(nat * 2)}, ${fmt(nat * 5)} or ${fmt(nat * 10)} m)`;
			}
		}
		return null;
	}
	function resCheck() {
		const e = resError();
		$('#s-res-err').textContent = e ? e[0].toUpperCase() + e.slice(1) + '.' : '';
		$('#s-res-err').classList.toggle('hidden', !e);
		return !e;
	}

	function values() {
		const cover = $('#s-cover').value;
		return {ECOSIM_CUT: $('#s-cut').value || '', ECOSIM_AREA: $('#s-area').value || '', ECOSIM_RESOLUTION: resConfig(),
			ECOSIM_STAT: $('#s-stat').value,
			ECOSIM_NEIGHBOURS: Number($('#s-nb').value), ECOSIM_MIN_COVER: cover === '' ? null : Number(cover),
			ECOSIM_ECOREGIONS: $('#s-eco').value.trim(), ECOSIM_MARINE: $('#s-marine').value.trim()};
	}

	function markChanged() {
		const v = values();
		let nd = 0, unsaved = 0;
		for (const p of S.params) {
			const d = String(v[p.name] ?? '') !== String(p.default);
			const f = $('#sf-' + p.name);
			if (f) f.classList.toggle('changed', d);
			if (d) nd++;
			if (S.saved && JSON.stringify(v[p.name]) !== JSON.stringify(S.saved[p.name])) unsaved++;
		}
		$('#s-changed').textContent = (unsaved ? `${unsaved} unsaved changes · ` : '') + (nd ? `${nd} options not at their default`
			: 'all options at their default');
	}

	// Starting HFTs ------------------------------------------------------------------------------------------------
	async function loadList() {
		try { S.list = await api('results_list'); } catch (e) { fail(e); S.list = []; }
		const sel = $('#s-sc');
		const cur = sel.value || store('ecosim.sc');
		sel.textContent = '';
		for (const x of S.list) {
			for (const sc of x.scenarios) {
				sel.append(el('option', {value: key(x.study_area, sc)},
					`${x.study_area} / ${sc}${sc.endsWith('_previous') ? ' (previous)' : ''}`));
			}
		}
		const has = v => [...sel.options].some(o => o.value === v);
		let want = cur;
		if (!want || !has(want)) {
			const st = await api('status').catch(() => null);
			want = st ? key(st.study_area, st.scenario) : null;
		}
		if (sel.options.length) sel.value = has(want) ? want : sel.options[0].value;
		const empty = !sel.options.length;
		$('#s-empty').classList.toggle('hidden', !empty);
		$('#s-src').classList.toggle('hidden', empty);
		return sel.value;
	}

	async function loadSources(selectRun) {
		const [sa, sc] = parts();
		if (!sa) {
			S.src = null; S.plan = null; S.planError = null;
			renderInfo(); renderPlan(); await loadRuns();
			return;
		}
		store('ecosim.sc', $('#s-sc').value);
		// the pixel size of the form (or the one waiting for the plan of the area chosen; at the start, that of config.py)
		const wantRes = S.resFollow !== null ? S.resFollow : S.src ? resConfig() : Number(S.cfg.ECOSIM_RESOLUTION) || 0;
		const wantCut = S.src ? $('#s-cut').value : null;
		const wantArea = S.src ? $('#s-area').value : S.cfg.ECOSIM_AREA || '';
		$('#s-plan-sub').textContent = 'loading…';
		// a plan still coming for the scenario chosen before must not replace the plan of this one
		clearTimeout(planTimer);
		planTimer = null;
		const seq = ++S.planSeq;
		let src;
		try {
			src = await api('ecosim_sources', sa, sc);
		} catch (e) {
			fail(e);
			src = null;
		}
		if (key(sa, sc) !== $('#s-sc').value) return; // another scenario was chosen meanwhile
		S.src = src; S.stale = false;
		if (src) {
			fillFields(src.options.params);
			S.cfg = Object.fromEntries(src.options.params.map(p => [p.name, p.value]));
			if (src.plan && src.plan.ecoregions && !$('#s-eco').value.trim()) {
				S.defaultEco = src.plan.ecoregions.file; S.defaultEcoPath = src.plan.ecoregions.path;
			}
			fillCuts(wantCut === null ? src.cut : wantCut);
			fillAreas(wantArea);
		}
		// the plan of ecosim_sources is made with the options of config.py, also its area and pixel size: it is the plan
		// of the form only when the form has that area; else (or when it failed) the area of the form is planned at its
		// own pixel, and the pixel size of the form is kept if it is a multiple of it (S.resFollow)
		const own = src && src.plan && $('#s-area').value === areaValue(S.cfg.ECOSIM_AREA) ? src.plan : null;
		const plannable = !!(src && src.ids.length && src.cuts.length);
		// a plan asked meanwhile (seq !== S.planSeq) was made with the form of before: it is asked again below
		let again = plannable && (!own || seq !== S.planSeq);
		S.plan = own; S.planError = null;
		if (own || !plannable) { S.resFollow = null; fillRes(keepRes(wantRes, own)); }
		else S.resFollow = wantRes ?? 0;
		ecoPlaceholder();
		markChanged();
		if (own && !again) {
			const v = values();
			again = !(v.ECOSIM_CUT === src.cut && String(v.ECOSIM_RESOLUTION) === String(Number(S.cfg.ECOSIM_RESOLUTION) || 0)
				&& String(v.ECOSIM_MIN_COVER) === String(S.cfg.ECOSIM_MIN_COVER) && v.ECOSIM_ECOREGIONS === (S.cfg.ECOSIM_ECOREGIONS || ''));
		}
		// without the plan of this area the page keeps what it shows until refreshPlan draws the new plan
		if (own || !again) { renderInfo(); renderPlan(); }
		const runs = loadRuns(selectRun);
		if (again) await refreshPlan();
		await runs;
	}

	function fillCuts(want) {
		const sel = $('#s-cut');
		sel.textContent = '';
		const src = S.src;
		if (!src) return;
		const multi = src.ids.length > 1;
		const sug = src.suggested || {};
		const sugCuts = [...new Set(Object.values(sug))];
		// one suggested cut for all the study areas: it is the value '' (ECOSIM_CUT = the suggested one)
		const one = sugCuts.length === 1 && Object.keys(sug).length === src.ids.length ? sugCuts[0] : null;
		const n = c => {
			const v = Object.values(c.hfts || {});
			return v.length ? ` · ${multi ? v.join(' / ') : v[0]} HFTs` : '';
		};
		if (sugCuts.length && !one) {
			sel.append(el('option', {value: ''},
				'The suggested cut of each study area (' + src.ids.map(i => `${i}: ${sug[i] || '—'}`).join(', ') + ')'));
		}
		for (const c of src.cuts) {
			const name = c.cut === 'results' ? 'Scenario folder (SIMIL 1 or 2)' : c.cut;
			sel.append(el('option', {value: c.cut === one ? '' : c.cut}, name + n(c) + (c.cut === one ? ' · suggested' : '')));
		}
		const has = v => [...sel.options].some(o => o.value === v);
		const w = want === one ? '' : want;
		sel.value = has(w) ? w : has(src.cut) ? src.cut : (sel.options[0] || {}).value || '';
	}

	function ecoPlaceholder() {
		$('#s-eco').placeholder = S.defaultEco ? `default: ${S.defaultEco}` : 'default: ecoregion/ecoregions.gpkg of the set of variables';
	}

	// the error of ecosim_sources (no HFTs, or the plan with the options of config.py): not while the plan of another
	// area is coming (S.resFollow)
	const srcError = () => S.src && !S.plan && S.resFollow === null ? S.src.error : null;

	function renderInfo() {
		const p = S.plan; // that of the area and options of the form (not the plan of config.py of ecosim_sources)
		const facts = $('#s-facts');
		facts.textContent = '';
		const fact = (k, v, s) => el('div', {class: 'fact'}, el('div', {class: 'k'}, k), el('div', {class: 'v'}, v),
			s ? el('div', {class: 's'}, s) : null);
		if (S.src) {
			const ids = S.src.ids;
			facts.append(fact('Study areas', `${ids.length} with HFTs`,
				ids.slice(0, 4).map(areaName).join(', ') + (ids.length > 4 ? ', …' : '')));
		}
		const ar = p && p.area;
		const comp = ar && ar.kind === 'comparison area';
		if (p && p.native) {
			const nt = p.native;
			facts.append(fact(comp ? 'Comparison area grid' : 'Ecoregion grid', `${fmt(nt.res)} m · ${fmt(nt.width, 0)} × ${fmt(nt.height, 0)} px`,
				String(nt.crs || '').replace(/ \(.*\)$/, '')));
		}
		if (ar) {
			const own = ar.reference && ar.reference !== ar.folder;
			facts.append(fact('HFT statistics from', own ? 'the study area at the pixel of the area' : 'the grid of the area',
				own ? ar.reference : 'its pixels covered by the HFTs'));
		}
		// what the HFTs are compared with, and from which grid their statistics come
		const ai = $('#s-area-info');
		ai.textContent = '';
		if (ar) {
			const outside = ((p && p.targets) || []).filter(t => t.in_area === false).map(t => t.id);
			const ref = ar.reference && ar.reference !== ar.folder;
			ai.append(el('b', {}, comp ? 'Comparison area' : 'Ecoregions'), `: ${ar.name || ar.folder}`,
				ar.tag ? el('span', {class: 'faint'}, ` (results in similarity_<pixel>m_${ar.tag})`) : '', '. ',
				ref ? `The statistics of the HFTs come from ${ar.reference}, the study area at the pixel and CRS of the area.`
					: 'The statistics of the HFTs come from the pixels of the area that they cover.',
				outside.length ? ` The study area${outside.length > 1 && S.src && S.src.ids.length > 1 ? 's ' + outside.join(', ') + ' are'
					: ' is'} outside the area: its HFTs are only compared with it.` : '');
		}
		const vars = (p && p.variables) || [];
		const paired = vars.filter(v => v.ecoregion);
		const missing = ((p && p.missing) || []).filter(k => !vars.some(v => v.key === k));
		if (p) {
			facts.append(fact('Variables', `${paired.length} of ${vars.length + missing.length} paired`,
				paired.length === vars.length && !missing.length ? 'plot scale ↔ ecoregion scale'
					: 'the others have no file of the same name in ecoregion'));
		}
		const tbl = $('#s-vars');
		tbl.textContent = '';
		tbl.closest('.table-wrap').classList.toggle('hidden', !p);
		if (p) {
			tbl.append(el('colgroup', {}, ...['32%', '26%', '26%', '16%'].map(w => el('col', {style: {width: w}}))));
			tbl.append(el('thead', {}, el('tr', {}, ...['Variable', 'Plot scale (study area)', comp ? 'Comparison area' : 'Ecoregion scale', '']
				.map(h => el('th', {}, h)))));
			const body = el('tbody');
			for (const v of vars) {
				body.append(el('tr', {}, el('td', {}, el('div', {class: 'vname'}, v.key)),
					el('td', {class: 'small', title: v.local || ''}, base(v.local)),
					el('td', {class: 'small', title: v.ecoregion || ''}, v.ecoregion ? base(v.ecoregion) : 'no file of this name'),
					el('td', {}, v.ecoregion ? el('span', {class: 'badge ok'}, 'OK') : el('span', {class: 'badge err'}, 'missing'))));
			}
			for (const k of missing) {
				body.append(el('tr', {}, el('td', {}, el('div', {class: 'vname'}, k)),
					el('td', {class: 'small', colspan: 2}, 'no file of this name at the other scale'),
					el('td', {}, el('span', {class: 'badge err'}, 'missing'))));
			}
			if (!vars.length && !missing.length) body.append(el('tr', {}, el('td', {colspan: 4, class: 'empty'}, 'No variables.')));
			tbl.append(body);
		}
		const info = $('#s-eco-info');
		info.textContent = '';
		const e = p && p.ecoregions;
		if (e) {
			const fs = e.features || [];
			info.append(`${fs.length} ${fs.length === 1 ? 'ecoregion' : 'ecoregions'} in ${base(e.file)}`
				+ `${e.field ? ' (column ' + e.field + ')' : ''}: `,
				fs.slice(0, 5).map(f => `${f.id}${f.name ? ' ' + f.name : ''}`).join('; ') + (fs.length > 5 ? '; …' : ''));
		} else if (p && p.native && !(p.errors || []).length) {
			info.append(el('span', {class: 'err-text'}, 'No polygons of the ecoregions: each HFT is compared with all the pixels of the '
				+ 'ecoregion variables.'));
		}
		const notes = $('#s-notes');
		notes.textContent = '';
		const err = S.planError || srcError();
		if (err) notes.append(el('div', {class: 'notice err'}, el('div', {}, el('b', {}, 'The similarity cannot be planned: '), err)));
		// the study_area folder of the area that similarity.py does not use (ecosim: plan.area.reference_note): said in the
		// notice of the HFTs without pixels when there is one, not twice
		const rejected = p && p.reference_rejected;
		const told = t => p.no_hft_pixels && rejected && String(t).endsWith('is not used: ' + rejected);
		for (const t of (p && p.errors) || []) notes.append(el('div', {class: 'notice err'}, t));
		for (const t of (p && p.warnings) || []) if (!told(t)) notes.append(el('div', {class: 'notice warn'}, t));
		if (p && p.no_hft_pixels) {
			const where = 'in "Download from MonitorEO" (with the same study area: the download adds ecoregion/study_area when '
				+ 'the area does not contain it), or choose another area.';
			const why = rejected
				? 'It does not cover the study area, and its study_area folder (the study area at its pixel) cannot be used: '
					+ `${String(rejected).replace(/\.\s*$/, '')}. Download the area again with the variables of this set ${where}`
				: 'It does not cover the study area and it has no study_area folder (the study area at its pixel). Download it '
					+ `again ${where}`;
			notes.append(el('div', {class: 'notice warn'}, el('div', {}, el('b', {}, 'The HFTs get no pixels from this area. '), why)));
		}
	}

	// HFT pixels at this pixel size -------------------------------------------------------------------------------
	let planTimer = null;
	function schedulePlan() {
		clearTimeout(planTimer);
		markChanged();
		resCheck();
		$('#s-plan-sub').textContent = 'updating…';
		planTimer = setTimeout(() => refreshPlan(), 450);
		canRun();
	}

	async function refreshPlan() {
		clearTimeout(planTimer);
		planTimer = null;
		const [sa, sc] = parts();
		if (!S.src || !sa || !S.src.ids.length) return;
		// another area (S.resFollow): planned at its own pixel, whatever the select has (the pixels of the area of before)
		const follow = S.resFollow !== null;
		if (!follow && !resCheck()) { $('#s-plan-sub').textContent = ''; canRun(); return; }
		const seq = ++S.planSeq;
		const cover = $('#s-cover').value === '' ? null : Number($('#s-cover').value);
		const eco = $('#s-eco').value.trim();
		const scKey = key(sa, sc);
		$('#s-plan-sub').textContent = 'updating…';
		try {
			// 0 = the pixel of the ecoregion variables (not "that of config.py")
			const p = await api('ecosim_plan', sa, sc, $('#s-cut').value || '', follow ? 0 : resMetres() ?? 0, cover, eco,
				$('#s-area').value || '');
			if (seq !== S.planSeq || scKey !== $('#s-sc').value) return;
			S.plan = p; S.planError = null;
			if (!eco && p.ecoregions) { S.defaultEco = p.ecoregions.file; S.defaultEcoPath = p.ecoregions.path; ecoPlaceholder(); }
			if (follow && planNative(p)) {
				// the select of this area, with the pixel size of before if it is a multiple of its pixel, else its pixel
				const want = keepRes(S.resFollow, p);
				S.resFollow = null;
				fillRes(want);
				markChanged();
				if (want) return await refreshPlan(); // the plan at the pixel size kept
			} else if (follow || resKey(p) !== S.resKey) {
				// a plan without the pixel of the area (it lacks variables: S.resFollow waits for one with it), or with
				// other pixel sizes than those of the select (e.g. the files of the area changed): the select of this plan
				fillRes(follow ? 0 : planNative(p) && Math.abs(p.res - planNative(p)) > 1e-9 ? p.res : 0);
				markChanged();
			}
		} catch (e) {
			if (seq !== S.planSeq || scKey !== $('#s-sc').value) return;
			S.plan = null; S.planError = e.message;
			if (follow) { fillRes(0); markChanged(); } // no pixel sizes of the area of before; S.resFollow waits
		}
		renderInfo(); renderPlan();
	}

	function renderPlan() {
		const box = $('#s-plan');
		box.textContent = '';
		const p = S.plan;
		const note = $('#s-plan-note');
		if (!p) {
			box.append(el('div', {class: 'empty'}, S.planError || srcError() || (S.src ? 'No plan.' : 'Choose the HFTs.')));
			$('#s-plan-sub').textContent = '';
			note.classList.add('hidden');
			canRun();
			return;
		}
		$('#s-plan-sub').textContent = !p.grid ? '' : `${fmt(p.res)} m pixels · ${fmt(p.grid.width, 0)} × ${fmt(p.grid.height, 0)} px`
			+ (p.memory_mb ? ` · ${mem(p.memory_mb)} of memory` : '');
		const multi = (p.targets || []).length > 1;
		for (const t of p.targets || []) {
			const hfts = t.hfts || [];
			const count = st => hfts.filter(h => h.status === st).length;
			const summary = [`${count('ok')} ok`];
			if (count('few_pixels')) summary.push(`${count('few_pixels')} with few pixels`);
			if (count('too_few_pixels')) summary.push(`${count('too_few_pixels')} without similarity (too few pixels)`);
			const tbl = el('table', {class: 't'});
			tbl.append(el('thead', {}, el('tr', {}, el('th', {}, 'HFT'), el('th', {class: 'num'}, 'Area (ha)'),
				el('th', {class: 'num', title: 'Pixels of the grid of the similarity that belong to the HFT'}, 'Pixels'),
				el('th', {}, 'Status'))));
			const body = el('tbody');
			for (const h of hfts) {
				body.append(el('tr', {class: h.status === 'too_few_pixels' ? 'off' : ''},
					el('td', {}, el('span', {class: 'hft'}, swatch(h.colour), h.hft)),
					el('td', {class: 'num'}, fmt(h.area_ha, 1)), el('td', {class: 'num'}, fmt(h.ref_pixels, 0)),
					el('td', {}, statusBadge(h.status))));
			}
			if (!hfts.length) body.append(el('tr', {}, el('td', {colspan: 4, class: 'empty'}, 'No HFTs.')));
			tbl.append(body);
			const cut = t.cut === 'results' ? 'scenario folder' : `cut ${t.cut}${t.suggested ? ' (suggested)' : ''}`;
			box.append(el('div', {class: 'stack', style: {gap: '6px'}},
				el('div', {class: 'plan-head'}, el('h3', {}, multi ? areaName(t.id) : cut[0].toUpperCase() + cut.slice(1)),
					el('span', {class: 'sub'}, (multi ? cut + ' · ' : '') + `${hfts.length} HFTs: ` + summary.join(', '))),
				el('div', {class: 'table-wrap', style: {maxHeight: '360px'}}, tbl)));
		}
		if (!(p.targets || []).length) {
			box.append(el('div', {class: 'empty'}, p.grid ? 'No HFTs for this cut.' : 'The pixels of the HFTs cannot be counted: see the '
				+ 'errors in "Starting HFTs".'));
		}
		const cover = Number($('#s-cover').value) || 0.5;
		note.textContent = '';
		note.classList.toggle('hidden', !p.grid);
		note.append(el('div', {class: 'stack', style: {gap: '6px'}},
			el('div', {}, el('b', {}, 'Why the grid of the ecoregion. '), 'The statistics of each HFT (the mean and covariance of its '
				+ 'variables) are taken from the ecoregion variables, on the pixels of this grid covered by the HFT, so that the HFT and '
				+ 'the rest of its ecoregion are compared with the same rasters and the same pixel size. The fine pixels of the HFT map '
				+ 'are not used: they vary more within an HFT than means over larger pixels, and they come from other rasters, so the '
				+ 'comparison would mix two scales.'),
			el('div', {}, el('b', {}, 'Small HFTs. '),
				`A pixel belongs to an HFT when the HFT covers at least ${fmt(cover * 100, 0)} % of it. `
				+ 'With coarser pixels a small or fragmented HFT covers fewer of them: '
				+ (p.min_pixels ? `with fewer than ${p.min_pixels} (the variables + 1, at least 5) its similarity is not computed`
					: 'with too few its similarity is not computed')
				+ (p.few_pixels ? `, and with fewer than ${p.few_pixels} (3 × the variables) it is less reliable` : '')
				+ '. A finer pixel size or a cut with fewer, larger HFTs keeps them.')));
		canRun();
	}

	// Run ---------------------------------------------------------------------------------------------------------
	function blockers() {
		if (view.running()) return [];
		const why = [];
		if (!S.src) why.push('choose the HFTs (a scenario with results)');
		else if (S.planError || srcError()) why.push(S.planError || srcError());
		for (const t of (S.plan && S.plan.errors) || []) why.push(t);
		if (S.plan && S.plan.grid && !(S.plan.targets || []).length) why.push('no HFTs for this cut');
		const re = resError();
		if (re) why.push(re);
		return why;
	}

	function canRun() {
		const why = blockers();
		$('#s-run').disabled = view.running() || why.length > 0 || !S.plan;
		const text = why.join('; ');
		$('#s-run-why').textContent = why.length ? 'Cannot run: ' + (text.length > 180 ? text.slice(0, 180) + '…' : text + '.') : '';
		$('#s-run-why').title = text;
		const sub = $('#s-run-sub');
		sub.textContent = '';
		const ts = (S.plan && S.plan.targets) || [];
		if (ts.length && !view.running()) {
			const old = ts.filter(t => t.run_exists).length;
			sub.append(`Results in ${ts[0].run_folder}${ts.length > 1 ? ` and ${ts.length - 1} more` : ''}`,
				old ? el('span', {class: 'badge warn', style: {marginLeft: '6px'}}, 'exists: it will be replaced') : '');
		}
	}

	async function pickVector(input, title, isDefault) {
		try {
			const f = await pickFiles('vector', false, title);
			if (!f || !f.length) return;
			input.value = isDefault && isDefault(f[0]) ? '' : f[0];
			schedulePlan();
		} catch (e) { fail(e); }
	}

	$('#s-sc').addEventListener('change', () => loadSources());
	$('#s-refresh').addEventListener('click', async () => { await loadList(); await loadSources(); });
	$('#s-cut').addEventListener('change', schedulePlan);
	// the pixel sizes of the select are those of the area of before: the new one is planned at its own pixel, and the
	// pixel size chosen is kept if it is a multiple of it (also one chosen before the plan of the new area arrives)
	$('#s-area').addEventListener('change', () => { if (S.resFollow === null) S.resFollow = resConfig() ?? 0; schedulePlan(); });
	const resChanged = () => { if (S.resFollow !== null) S.resFollow = resConfig() ?? 0; schedulePlan(); };
	$('#s-res').addEventListener('change', () => { syncRes(); resChanged(); });
	$('#s-res-other').addEventListener('input', resChanged);
	$('#s-cover').addEventListener('input', schedulePlan);
	$('#s-eco').addEventListener('change', schedulePlan);
	$('#s-stat').addEventListener('change', markChanged);
	$('#s-nb').addEventListener('change', markChanged);
	$('#s-marine').addEventListener('input', markChanged);
	$('#s-eco-pick').addEventListener('click', () => pickVector($('#s-eco'), 'File of the ecoregions', f => f === S.defaultEcoPath));
	$('#s-eco-clear').addEventListener('click', () => { $('#s-eco').value = ''; schedulePlan(); });
	$('#s-marine-pick').addEventListener('click', () => pickVector($('#s-marine'), 'File of the marine ecoregions'));
	$('#s-marine-clear').addEventListener('click', () => { $('#s-marine').value = ''; markChanged(); });
	$('#s-save').addEventListener('click', async () => {
		if (S.resFollow !== null) await refreshPlan(); // the pixel size of the area chosen, not of the area of before
		if (resError()) return toast(resError()[0].toUpperCase() + resError().slice(1) + '.', 'warn');
		try {
			const ch = await api('ecosim_save', values());
			toast(ch.length ? `Saved in config.py: ${ch.join(', ')}.` : 'There were no changes.', 'ok');
			S.saved = values();
			markChanged();
		} catch (e) { fail(e); }
	});
	$('#s-defaults').addEventListener('click', async () => {
		if (!await confirmBox('Default values', 'The form gets the default values of all the options (not saved yet).')) return;
		setValues(S.defaults);
		S.resFollow = Number(S.defaults.ECOSIM_RESOLUTION) || 0; // the area may be another one: its own pixel
		schedulePlan();
	});
	$('#s-run').addEventListener('click', async () => {
		if (view.running()) return;
		const [sa, sc] = parts();
		if (!sa) return;
		await refreshPlan();
		const why = blockers();
		if (why.length || !S.plan) return toast('The similarity cannot be computed: ' + (why.join('; ') || 'no plan') + '.', 'warn');
		const old = (S.plan.targets || []).filter(t => t.run_exists);
		if (old.length) {
			const ok = await confirmBox('The similarity already has results', el('div', {class: 'stack', style: {gap: '8px'}},
				el('div', {}, 'These folders already have results of the similarity with this pixel size and area:'),
				el('ul', {class: 'small'}, ...old.map(t => el('li', {}, t.run_folder))),
				el('div', {}, 'If you go on, they are replaced. To keep them, choose another pixel size or copy the folder first.')),
				'Replace');
			if (!ok) return;
		}
		try {
			const id = await api('ecosim_run', values(), sa, sc, $('#s-cut').value || '');
			S.saved = values();
			markChanged();
			view.start(id, TITLE);
		} catch (e) { fail(e); }
	});
	$('#s-stop').addEventListener('click', async () => {
		const text = 'The run stops; its folder is left unfinished (without similarity_run.json).';
		if (await confirmBox('Stop the similarity', text, 'Stop')) stopJob(view);
	});

	// Results -----------------------------------------------------------------------------------------------------
	function runLabel(r) {
		const c = r.counts || {};
		const without = r.hfts - (c.ok || 0) - (c.few_pixels || 0) - (c.no_similar || 0);
		// another area than the ecoregion of the set: its tag (the suffix of the folder) and name
		const a = r.area || {};
		const area = r.tag ? ` · vs ${a.name ? a.name + ' (' + r.tag + ')' : r.tag}` : '';
		return `${r.cut === 'results' ? 'Scenario folder' : r.cut} · ${fmt(r.res)} m${area} · ${r.hfts} HFTs`
			+ (without > 0 ? ` (${without} without similarity)` : '') + (r.date ? ' · ' + r.date : '');
	}

	async function loadRuns(select) {
		const [sa, sc] = parts();
		S.runs = sa ? await api('ecosim_runs', sa, sc).catch(e => { fail(e); return []; }) : [];
		S.previews.clear();
		const sel = $('#s-runs');
		const cur = select || sel.value;
		sel.textContent = '';
		for (const r of S.runs) sel.append(el('option', {value: r.folder}, runLabel(r)));
		const empty = !S.runs.length;
		$('#s-res-empty').classList.toggle('hidden', !empty);
		for (const id of ['#s-res-body', '#s-runs', '#s-open', '#s-map-card', '#s-plots-card']) $(id).classList.toggle('hidden', empty);
		$('label[for=s-runs]').classList.toggle('hidden', empty);
		if (empty) {
			S.res = null;
			$('#s-res-sub').textContent = '';
			return;
		}
		sel.value = S.runs.some(r => r.folder === cur) ? cur : S.runs[0].folder;
		await loadResult();
	}

	async function loadResult() {
		const folder = $('#s-runs').value;
		if (!folder) return;
		let R;
		try { R = await api('ecosim_results', folder); } catch (e) { return fail(e); }
		if (folder !== $('#s-runs').value) return;
		S.res = R;
		const run = R.run || {};
		const o = run.options || {};
		const ar = R.area || {};
		const vs = ar.tag ? 'compared with ' + (ar.kind === 'comparison area' ? 'the comparison area ' : '') + (ar.name || ar.tag) : null;
		$('#s-res-sub').textContent = [vs, run.res ? `${fmt(run.res)} m pixels` : null, o.stat ? `threshold: ${o.stat}` : null,
			o.neighbours ? `${o.neighbours} neighbours` : null, o.min_cover !== undefined ? `minimum cover ${o.min_cover}` : null,
			(run.variables || []).length ? `${run.variables.length} variables` : null,
			run.seconds ? (run.seconds >= 60 ? fmt(run.seconds / 60, 1) + ' min' : fmt(run.seconds, 0) + ' s') : null]
			.filter(Boolean).join(' · ');
		renderThdi(R);
		renderHri(R);
		$('#s-meta').textContent = R.metadata || 'No similarity_metadata.txt';
		gallery($('#s-plots'), R.plots);
		const multi = new Set(R.maps.map(m => m.id)).size > 1;
		const hsel = $('#s-map-hft');
		const cur = hsel.value;
		hsel.textContent = '';
		for (const m of R.maps) {
			hsel.append(el('option', {value: `${m.id}_${m.hft}`}, `${multi ? 'ID ' + m.id + ' · ' : ''}HFT ${m.hft}`
				+ (m.sim || m.similar ? '' : ` (${(STATUS[m.status] || [, m.status])[1]})`)));
		}
		const first = R.maps.find(m => m.sim || m.similar) || R.maps[0];
		hsel.value = R.maps.some(m => `${m.id}_${m.hft}` === cur) ? cur : first ? `${first.id}_${first.hft}` : '';
		await drawMap();
	}

	function renderThdi(R) {
		const box = $('#s-thdi');
		box.textContent = '';
		if (!R.thdi.length) {
			box.append(el('div', {class: 'small muted'}, 'No THDI of the study areas (thdi_res.csv).'));
			return;
		}
		for (const row of R.thdi) {
			if (R.thdi.length > 1) box.append(el('div', {class: 'thdi-id'}, areaName(row.wdpaid)));
			box.append(el('div', {class: 'facts'}, ...THDI.map(([k, label, s]) => el('div', {class: 'fact', title: k},
				el('div', {class: 'k'}, label), el('div', {class: 'v'}, fmt(row[k], 3)), s ? el('div', {class: 's'}, s) : null))));
		}
	}

	function renderHri(R) {
		const tbl = $('#s-hri');
		tbl.textContent = '';
		const t = R.tables.hri;
		$('#s-hri-file').textContent = t ? t.file : '';
		fullTable($('#s-hri-all'), t);
		if (!t) {
			tbl.append(el('tbody', {}, el('tr', {}, el('td', {class: 'empty'}, 'No hri_results.csv.'))));
			return;
		}
		const tcol = t.header.includes('MedianHFTSim') ? 'MedianHFTSim' : 'AveHFTSim';
		const multi = new Set(R.hri.map(r => r.wdpaid)).size > 1;
		const colour = Object.fromEntries(R.maps.map(m => [`${m.id}_${m.hft}`, m.colour]));
		const cols = HRI.filter(([c]) => c === 'T' || t.header.includes(c));
		tbl.append(el('thead', {}, el('tr', {}, multi ? el('th', {}, 'ID') : null, el('th', {}, 'HFT'), el('th', {}, 'Status'),
			...cols.map(([c, h, , title]) => el('th', {class: 'num', title: c === 'T' ? `${title} (${tcol})` : `${title} (${c})`}, h)))));
		const body = el('tbody');
		for (const r of R.hri) {
			const hft = r.hft ?? String(r.wdpaid_hft || '').split('_').pop();
			const k = `${r.wdpaid}_${hft}`;
			const tr = el('tr', {class: 'click', title: 'Show its maps'},
				multi ? el('td', {}, r.wdpaid) : null,
				el('td', {}, el('span', {class: 'hft'}, swatch(colour[k]), hft)),
				el('td', {}, statusBadge(r.Status)),
				...cols.map(([c, , d]) => el('td', {class: 'num'}, fmt(c === 'T' ? r[tcol] : r[c], d))));
			tr.addEventListener('click', () => {
				if (![...$('#s-map-hft').options].some(o => o.value === k)) return;
				$('#s-map-hft').value = k;
				drawMap();
				$('#s-map-card').scrollIntoView({behavior: 'smooth', block: 'start'});
			});
			body.append(tr);
		}
		tbl.append(body);
	}

	function fullTable(tbl, t) {
		tbl.textContent = '';
		if (!t) return;
		tbl.append(el('thead', {}, el('tr', {}, ...t.header.map(h => el('th', {}, h)))));
		const body = el('tbody');
		for (const r of t.rows) {
			body.append(el('tr', {}, ...r.map(v => el('td', {class: /^-?[\d.]+(e-?\d+)?$/.test(v) ? 'num' : ''},
				/^-?\d+\.\d{4,}/.test(v) ? fmt(v, 4) : v))));
		}
		tbl.append(body);
	}

	// Map ---------------------------------------------------------------------------------------------------------
	function ensureMap() {
		if (!S.map) {
			S.map = makeMap('s-map', 'similarity', 'osm', $('#s-base'));
			S.map.on('moveend', () => { clearTimeout(S.detailTimer); S.detailTimer = setTimeout(detail, 300); });
		}
		setTimeout(() => S.map.invalidateSize(), 60);
	}

	const opacity = () => Number($('#s-op').value) / 100;

	// the whole map is at most 1200 pixels wide: zoomed in beyond them, the pixels of the view are drawn over it
	async function detail() {
		const p = S.full;
		if (!S.map || !p || !S.overlay) return dropDetail();
		const v = S.map.getBounds();
		const [w0, s0, e0, n0] = p.bbox;
		const shown = (v.getEast() - v.getWest()) / ((e0 - w0) / p.width); // pixels of the whole map across the view
		if (!bboxBounds(p.bbox).intersects(v) || S.map.getSize().x / shown < 1.5) return dropDetail();
		const padX = (v.getEast() - v.getWest()) * 0.25, padY = (v.getNorth() - v.getSouth()) * 0.25;
		const box = [Math.max(w0, v.getWest() - padX), Math.max(s0, v.getSouth() - padY), Math.min(e0, v.getEast() + padX),
			Math.min(n0, v.getNorth() + padY)];
		const seq = ++S.detailSeq;
		const want = S.fullKey;
		let d;
		try { d = await api('ecosim_preview', S.fullTif, S.fullKind, box); } catch (e) { return; }
		if (seq !== S.detailSeq || want !== S.fullKey || !S.overlay) return;
		if (S.detail) S.map.removeLayer(S.detail);
		S.detail = L.imageOverlay(d.image, bboxBounds(d.bbox), {opacity: opacity(), className: 'pixelated'}).addTo(S.map);
		S.overlay.setOpacity(0);
		coarseNote();
		front();
	}
	function dropDetail() {
		S.detailSeq++;
		if (S.detail) { S.map.removeLayer(S.detail); S.detail = null; }
		if (S.overlay) S.overlay.setOpacity(opacity());
		coarseNote();
	}
	// the note of the legend on the coarser whole map, not while the pixels of the view are shown
	function coarseNote() { $$('#s-legend .note').forEach(n => n.classList.toggle('hidden', !!S.detail)); }
	function front() {
		if (S.eco && S.map.hasLayer(S.eco)) S.eco.bringToFront();
		if (S.area && S.map.hasLayer(S.area)) S.area.bringToFront();
		if (S.saLayer && S.map.hasLayer(S.saLayer)) S.saLayer.bringToFront();
	}

	async function drawMap() {
		const R = S.res;
		if (!R) return;
		ensureMap();
		const m = R.maps.find(x => `${x.id}_${x.hft}` === $('#s-map-hft').value);
		const kind = $('#s-map-kind').value;
		const tif = m ? m[kind] : null;
		const multi = new Set(R.maps.map(x => x.id)).size > 1;
		$('#s-map-title').textContent = !m ? 'Similarity map'
			: `${kind === 'sim' ? 'Similarity' : 'Similar areas'} of HFT ${m.hft}${multi ? ' · ' + areaName(m.id) : ''}`;
		const seq = ++S.mapSeq;
		let p = null;
		if (tif) {
			$('#s-map-sub').textContent = 'loading…';
			const k = tif + '|' + kind;
			try {
				p = S.previews.get(k);
				if (!p) {
					p = await api('ecosim_preview', tif, kind);
					S.previews.set(k, p);
					if (S.previews.size > 16) S.previews.delete(S.previews.keys().next().value);
				}
			} catch (e) {
				if (seq === S.mapSeq) $('#s-map-sub').textContent = '';
				return fail(e);
			}
			if (seq !== S.mapSeq) return;
		}
		dropDetail();
		if (S.overlay) S.map.removeLayer(S.overlay);
		S.overlay = p ? L.imageOverlay(p.image, bboxBounds(p.bbox), {opacity: opacity(), className: 'pixelated'}).addTo(S.map) : null;
		S.full = p;
		S.fullKey = p ? tif + '|' + kind : null;
		S.fullTif = tif; S.fullKind = kind;
		$('#s-map-sub').textContent = m ? [m.threshold !== null && m.threshold !== undefined ? `T = ${fmt(m.threshold, 3)}` : null,
			m.ref_pixels !== null && m.ref_pixels !== undefined ? `${fmt(m.ref_pixels, 0)} HFT pixels` : null,
			m.area_ha ? `${fmt(m.area_ha, 1)} ha` : null, m.ecoregion ? `ecoregion ${m.ecoregion}` : null].filter(Boolean).join(' · ') : '';
		// the outlines of a run: again when the run was made again in the same folder (with other polygons, maybe)
		const ek = outlineKey(R);
		if (S.ecoFor !== ek) await ecoLayer(R, ek);
		legend(m, kind, p);
		await saLayer(S.src ? S.src.study_area : parts()[0]);
		front();
		// fitted only when the map has a size (not while the page is hidden: a run that ends on another page)
		if (S.fitted !== R.folder) { if (fitEco()) S.fitted = R.folder; } else detail();
	}
	function outlineKey(R) {
		return [R.folder, (R.ecoregions && R.ecoregions.path) || '', (R.area && R.area.path) || '', (R.run && R.run.date) || ''].join('|');
	}

	function legend(m, kind, p) {
		const box = $('#s-legend');
		box.textContent = '';
		box.classList.toggle('hidden', !m);
		if (!m) return;
		box.append(el('h4', {class: 'row', style: {gap: '6px'}}, swatch(m.colour), `HFT ${m.hft}`, statusBadge(m.status)));
		if (!p) {
			const why = m.status && m.status !== 'ok' ? ': ' + (STATUS[m.status] || [, , m.status])[2] : ' of this HFT';
			box.append(el('div', {class: 'small muted', style: {maxWidth: '230px'}}, `No map${why}.`));
		} else if (kind === 'sim') {
			const tick = m.threshold !== null && m.threshold !== undefined
				? el('div', {title: 'Threshold T of the HFT', style: {position: 'absolute', top: '-3px', bottom: '-3px', width: '2px',
					background: 'var(--ink)', left: `calc(${Math.max(0, Math.min(1, m.threshold)) * 100}% - 1px)`}}) : null;
			box.append(el('div', {class: 'small muted'}, 'Similarity (χ² probability)'),
				el('div', {style: {position: 'relative', width: '210px', margin: '4px 0 2px'}},
					el('div', {class: 'ramp', style: {background: `linear-gradient(90deg, ${p.ramp.join(',')})`}}), tick),
				el('div', {class: 'ramp-labels', style: {width: '210px'}}, el('span', {}, '0'), el('span', {}, '0.5'), el('span', {}, '1')),
				el('div', {class: 'legend-cat'}, el('span', {class: 'sw', style: {background: p.zero}}), el('span', {}, '0 (p ≤ 0.001)')),
				tick ? el('div', {class: 'small muted'}, `Bar: threshold T = ${fmt(m.threshold, 3)}`) : null);
		} else {
			for (const l of p.legend) {
				box.append(el('div', {class: 'legend-cat'}, el('span', {class: 'sw', style: {background: l.colour}}),
					el('span', {}, l.label)));
			}
		}
		if (p && p.note) {
			box.append(el('div', {class: 'small faint note' + (S.detail ? ' hidden' : ''), style: {maxWidth: '230px'}},
				'Coarser view: ' + p.note + '.'));
		}
		const line = (border, label) => el('div', {class: 'legend-cat'}, el('span', {class: 'sw', style: {background: 'transparent', border}}),
			el('span', {}, label));
		box.append(line('2px solid #ff3b30', 'Study area'));
		if (S.area) box.append(line(`2px dashed ${AREA_COLOUR}`, 'Comparison area'));
		if (S.eco) box.append(line('2px dashed #c99700', 'Ecoregion'));
	}

	const AREA_COLOUR = '#c2187a';
	const tip = (f, layer) => {
		const q = f.properties || {};
		if (q.id === undefined && !q.label) return;
		layer.bindTooltip(leafletText(`${q.id ?? ''}${q.label ? ' · ' + q.label : ''}`), {sticky: true});
	};
	// the ecoregions of the run and, when the HFTs were compared with a comparison area, its polygons
	async function ecoLayer(R, ek) {
		for (const k of ['eco', 'area']) if (S[k]) { S.map.removeLayer(S[k]); S[k] = null; }
		S.ecoFor = ek;
		const e = R.ecoregions;
		const a = R.area && R.area.kind === 'comparison area' && R.area.path && R.area.path !== (e && e.path) ? R.area : null;
		$('#s-area-show-lab').classList.toggle('hidden', !a);
		$('#s-fit-eco').textContent = a ? 'Zoom to the comparison area' : 'Zoom to the ecoregion';
		if (e && e.path) {
			const ids = (e.features || []).map(f => f.id).filter(x => x !== null && x !== undefined);
			const gj = await api('ecosim_outline', e.path, e.field || null, ids, R.folder).catch(err => { fail(err); return null; });
			if (S.ecoFor !== ek) return;
			if (gj) {
				S.eco = L.geoJSON(gj, {style: {color: '#c99700', weight: 2, dashArray: '6 4', fill: false}, onEachFeature: tip});
				if ($('#s-eco-show').checked) S.eco.addTo(S.map);
			}
		}
		if (a) {
			const gj = await api('ecosim_outline', a.path, null, null, R.folder).catch(err => { fail(err); return null; });
			if (S.ecoFor !== ek) return;
			if (gj) {
				S.area = L.geoJSON(gj, {style: {color: AREA_COLOUR, weight: 2.2, dashArray: '6 4', fill: false}, onEachFeature: tip});
				if ($('#s-area-show').checked) S.area.addTo(S.map);
			}
		}
	}

	async function saLayer(sa) {
		if (S.saFor === sa) return;
		if (S.saLayer) { S.map.removeLayer(S.saLayer); S.saLayer = null; }
		S.saFor = sa;
		const gj = sa ? await api('study_outline', sa).catch(() => null) : null;
		if (!gj || S.saFor !== sa) return;
		S.saLayer = L.geoJSON(gj, {style: {color: '#ff3b30', weight: 2.2, fill: false}, interactive: false});
		if ($('#s-sa-show').checked) S.saLayer.addTo(S.map);
	}

	function fitEco() { // to the comparison area, else the ecoregion, else the map; false when the map has no size (hidden)
		if (!S.map) return false;
		S.map.invalidateSize();
		const z = S.map.getSize();
		if (!z.x || !z.y) return false;
		const b = S.area ? S.area.getBounds() : S.eco ? S.eco.getBounds() : S.overlay ? S.overlay.getBounds() : null;
		if (b && b.isValid()) S.map.fitBounds(b, {padding: [12, 12]});
		return true;
	}
	function fitSa() {
		if (!S.map) return;
		S.map.invalidateSize();
		if (S.saLayer) S.map.fitBounds(S.saLayer.getBounds(), {padding: [30, 30], maxZoom: 14});
	}

	$('#s-runs').addEventListener('change', loadResult);
	$('#s-open').addEventListener('click', () => S.res && api('open_folder', S.res.path).catch(fail));
	$('#s-map-hft').addEventListener('change', drawMap);
	$('#s-map-kind').addEventListener('change', drawMap);
	$('#s-op').addEventListener('input', () => { (S.detail || S.overlay) && (S.detail || S.overlay).setOpacity(opacity()); });
	const toggle = (box, layer) => { if (layer) box.checked ? layer.addTo(S.map) : S.map.removeLayer(layer); };
	$('#s-eco-show').addEventListener('change', () => { toggle($('#s-eco-show'), S.eco); front(); });
	$('#s-area-show').addEventListener('change', () => { toggle($('#s-area-show'), S.area); front(); });
	$('#s-sa-show').addEventListener('change', () => { toggle($('#s-sa-show'), S.saLayer); front(); });
	$('#s-fit-sa').addEventListener('click', fitSa);
	$('#s-fit-eco').addEventListener('click', fitEco);

	async function init() {
		try {
			const o = await api('ecosim_options');
			fillFields(o.params);
			S.cfg = Object.fromEntries(o.params.map(p => [p.name, p.value]));
			setValues(S.cfg);
		} catch (e) { fail(e); }
		await loadList();
		await loadSources();
		S.saved = values();
		markChanged();
	}

	return {
		async onShow() {
			if (!S.loaded) {
				S.loaded = true;
				await init();
			} else {
				const before = $('#s-sc').value;
				await loadList();
				if ($('#s-sc').value !== before || S.stale) await loadSources();
				else await loadRuns();
			}
			if (S.map) setTimeout(() => S.map.invalidateSize(), 60);
		},
		attach(job) { view.start(job.id, job.title); },
		// the HFTs or the areas changed (a model run, a download): read them again now if the page is shown, else when
		// it is shown
		reload() {
			S.stale = true;
			if (S.loaded && $('#page-similarity').classList.contains('active')) this.onShow();
		},
	};
})();

// Updates from GitHub ================================================================================================
const updates = (() => {
	const view = new JobView(null, null, {onEnd: () => refresh()});
	const banner = $('#upd-banner');
	let dismissed = false;

	function show(kind, text, buttons = []) {
		banner.className = 'update-banner ' + kind;
		banner.textContent = '';
		banner.append(el('span', {}, text), el('div', {class: 'btns'}, ...buttons,
			el('button', {class: 'btn small ghost', onclick: () => { dismissed = true; banner.classList.add('hidden'); }}, '✕')));
	}

	async function refresh() {
		let st;
		try { st = await api('update_status'); } catch (e) { return; }
		const c = st.commit ? st.commit.slice(0, 7) : 'unknown';
		$('#upd-version').textContent = `${c} · ${st.branch || ''}`;
		$('#upd-version').title = `${st.repo} (${st.mode === 'git' ? 'git clone' : 'copy without git'})`;
		$('#upd-auto').checked = !!st.auto;
		const envRunning = st.env_job && (await api('running').catch(() => [])).some(j => j.id === st.env_job);
		if (envRunning) {
			show('warn', 'Updating the environment of eHabitat+ (new packages of environment.yml)…');
			setTimeout(refresh, 5000);
		} else if (st.restart) {
			show('ok', `eHabitat+ ${st.result && st.result.message ? st.result.message : 'was updated'}. Restart the app to use the new version.`,
				[el('button', {class: 'btn small primary', onclick: restart}, 'Restart now')]);
		} else if (view.shown && st.error) {
			show('warn', `Could not check for updates: ${st.error}`);
			view.shown = false;
		} else if (view.shown && st.result) {
			show('', `eHabitat+: ${st.result.message}.`);
			view.shown = false;
		} else if (st.just_updated && !dismissed) {
			show('ok', `eHabitat+ was updated when it opened: ${st.just_updated}.`);
		}
	}

	async function restart() {
		try {
			await api('restart');
			banner.textContent = 'Opening the new version…';
		} catch (e) { fail(e); }
	}

	$('#upd-check').addEventListener('click', async () => {
		if (view.running()) return;
		try {
			dismissed = false;
			show('', 'Checking for updates in GitHub…');
			view.shown = true;
			view.start(await api('check_updates'), 'Update');
		} catch (e) { fail(e); }
	});
	$('#upd-auto').addEventListener('change', e => api('set_auto_update', e.target.checked).catch(fail));
	return {refresh};
})();

// Start ==============================================================================================================
(async function boot() {
	try {
		const info = await api('info');
		ROOT = info.root;
		$('#root-path').textContent = info.root;
		$('#root-path').title = 'Working folder (copy of eHabpy)';
	} catch (e) { fail(e); }
	api('tools').then(t => {
		const box = $('#tools');
		box.textContent = '';
		for (const [k, label] of [['grass', 'GRASS'], ['r', 'R'], ['gdal', 'GDAL'], ['monitoreo', 'MonitorEO'], ['python', 'Python']]) {
			const v = t[k];
			box.append(el('span', {class: 'k'}, label), el('span', {class: 'v' + (v ? '' : ' bad'), title: v || 'not found'}, v ? String(v).replace(/^GRASS\s*/, '') : 'not found'));
		}
	}).catch(() => {});
	updates.refresh();
	await Promise.all([pages.inputs.init().catch(fail), pages.model.init().catch(fail), pages.gee.init().catch(fail)]);
	// jobs that were running when the interface was reloaded
	try {
		for (const j of await api('running')) {
			if (j.kind === 'model') pages.model.attach(j);
			else if (j.kind.startsWith('gee')) pages.gee.attach(j);
			else if (j.kind === 'ecosim') pages.similarity.attach(j);
		}
	} catch (e) { /* nothing */ }
	// MonitorEO at its last version of GitHub, in the background, when the app opens (if the app updates itself)
	api('gee_auto_update').then(id => {
		if (id) pages.gee.attach({id, kind: 'gee-update', title: 'Automatic update of MonitorEO'});
	}).catch(() => {});
	const p = store('page');
	if (p && pages[p] && p !== 'inputs') show(p);
})();
