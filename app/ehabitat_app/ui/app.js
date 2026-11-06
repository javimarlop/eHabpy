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
		else if (k === 'html') e.innerHTML = v;
		else if (v === true) e.setAttribute(k, '');
		else e.setAttribute(k, v);
	}
	for (const c of children.flat()) {
		if (c === null || c === undefined || c === false) continue;
		e.append(c instanceof Node ? c : document.createTextNode(String(c)));
	}
	return e;
}

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
				this.step = (/python|Rscript/i.test(w[0]) && w[1] ? w[1] : w[0]).split(/[\\/]/).pop();
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
	}});

	function ensureMap() {
		if (!S.map) S.map = makeMap('in-map', 'inputs', 'osm');
		setTimeout(() => S.map.invalidateSize(), 50);
	}

	// Study area ----------------------------------------------------------------------------------------------
	function renderDraft(info, ids, idError) {
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
		const facts = $('#sa-facts');
		facts.textContent = '';
		if (info) {
			facts.append(
				el('div', {class: 'fact'}, el('div', {class: 'k'}, 'Polygons'), el('div', {class: 'v'}, fmt(info.count, 0))),
				el('div', {class: 'fact'}, el('div', {class: 'k'}, 'Area'), el('div', {class: 'v'}, info.area_km2 !== null ? fmt(info.area_km2) + ' km²' : '—')),
				el('div', {class: 'fact'}, el('div', {class: 'k'}, 'CRS'), el('div', {class: 'v small'}, info.crs || 'no CRS')));
		}
		const chips = $('#sa-ids');
		chips.textContent = '';
		const sel = new Set(d && d.ids ? d.ids : ids);
		if (ids.length > 300) {
			chips.append(el('span', {class: 'small muted'}, `${ids.length} IDs: all of them are processed.`));
		} else {
			for (const i of ids) {
				const c = el('span', {class: 'chip' + (sel.has(i) ? ' on' : ''), 'data-id': i}, String(i));
				c.addEventListener('click', () => { c.classList.toggle('on'); dirty(); });
				chips.append(c);
			}
		}
		$('#sa-error').classList.toggle('hidden', !idError);
		$('#sa-error').textContent = idError || '';
	}

	function draftIds() {
		const on = $$('#sa-ids .chip.on').map(c => Number(c.dataset.id));
		return $$('#sa-ids .chip').length ? on : null;
	}

	function dirty() {
		const d = S.draft;
		const saved = S.sa;
		const same = saved && d && saved.file === d.file && (saved.layer || null) === (d.layer || null) && saved.col_id === d.col_id
			&& JSON.stringify((saved.ids || []).slice().sort()) === JSON.stringify((draftIds() || saved.ids || []).slice().sort());
		$('#sa-apply').disabled = !d || same;
		$('#sa-state').textContent = !d ? '' : same ? 'Used by the model.' : 'Changes not applied yet.';
	}

	async function loadStudyArea() {
		try {
			S.sa = await api('study_area');
		} catch (e) { fail(e); S.sa = null; }
		if (S.sa) {
			S.draft = {file: S.sa.file, layer: S.sa.info.layer, col_id: S.sa.col_id || 'cat', ids: S.sa.ids};
			S.sa.ids = S.sa.ids || S.sa.all_ids;
			renderDraft(S.sa.info, S.sa.all_ids, S.sa.info.id_error);
		} else {
			S.draft = null;
			renderDraft(null, [], null);
		}
		dirty();
		drawStudyArea();
	}

	async function inspect(file, layer, col) {
		const r = await api('inspect_vector', file, layer || null, col || 'cat');
		S.draft = {file, layer: r.info.layer, col_id: col || 'cat', ids: r.ids};
		renderDraft(r.info, r.ids, r.id_error);
		dirty();
		drawGeojson(await api('vector_preview', file, r.info.layer || null).catch(() => null));
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
		S.saLayer = L.geoJSON(gj, {style: {color: '#ff3b30', weight: 2.5, fill: true, fillOpacity: 0.06}}).addTo(S.map);
		S.map.fitBounds(S.saLayer.getBounds(), {padding: [20, 20]});
	}

	$('#sa-pick').addEventListener('click', async () => {
		try {
			const f = await pickFiles('vector', false, 'File of the study area');
			if (f && f.length) await inspect(f[0], null, 'cat');
		} catch (e) { fail(e); }
	});
	$('#sa-layer').addEventListener('change', () => inspect(S.draft.file, $('#sa-layer').value, 'cat').catch(fail));
	$('#sa-col').addEventListener('change', () => inspect(S.draft.file, S.draft.layer, $('#sa-col').value).catch(fail));
	$('#sa-all').addEventListener('click', () => { $$('#sa-ids .chip').forEach(c => c.classList.add('on')); dirty(); });
	$('#sa-none').addEventListener('click', () => { $$('#sa-ids .chip').forEach(c => c.classList.remove('on')); dirty(); });
	$('#sa-apply').addEventListener('click', async () => {
		const d = S.draft;
		const ids = draftIds();
		if (ids && !ids.length) return toast('Choose at least one ID.', 'warn');
		try {
			await api('set_study_area', d.file, d.layer || null, d.col_id, ids);
			toast('Study area applied (copied to segm_optim/pa, palist.csv written).', 'ok');
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
			s.append(el('option', {value: x.folder}, `${x.label} · ${x.study} + ${x.ecoregion} eco${x.folder === r.active ? ' · in use' : ''}`));
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
		if (S.folder === S.active) {
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

	function table(tbl, list, scale) {
		tbl.textContent = '';
		tbl.append(el('thead', {}, el('tr', {}, el('th', {}, ''), el('th', {}, 'Variable'), el('th', {class: 'num'}, 'Pixel'),
			el('th', {}, 'Size'), el('th', {}, 'NoData'), el('th', {}, 'CRS'), el('th', {}, ''))));
		const body = el('tbody');
		if (!list.length) body.append(el('tr', {}, el('td', {colspan: 7, class: 'empty'}, 'No variables.')));
		for (const r of list) {
			const cb = el('input', {type: 'checkbox', title: r.included ? 'Leave out of the set' : 'Include again'});
			cb.checked = r.included;
			cb.addEventListener('change', async () => {
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
		renderActive();
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
			S.eco = L.geoJSON(gj, {style: {color: '#ffcc00', weight: 2, dashArray: '6 4', fill: false}}).addTo(S.map);
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
	const S = {recipe: null, catalog: null};
	const view = new JobView($('#g-log'), $('#g-bar'), {
		onStart: () => { $('#g-stop').classList.remove('hidden'); $('#dot-gee').classList.remove('hidden'); },
		onEnd: s => {
			$('#g-stop').classList.add('hidden'); $('#dot-gee').classList.add('hidden');
			if (s.state === 'done' && view.kind === 'download') {
				toast('Download finished.', 'ok');
				pages.inputs.refresh();
				pages.model.refresh();
			}
			if (view.kind === 'update' || view.kind === 'connect' || view.kind === 'auth') status();
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
		const wc = new Set((rec.exclude_landcover || []).map(Number));
		for (const c of $$('#r-wc .chip')) c.classList.toggle('on', wc.has(Number(c.dataset.v)));
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
		rec.drive_folder = $('#r-drive-folder').value.trim() || 'eHabitat_GEE';
		rec.drive_dir = $('#r-drive-dir').value.trim() || null;
		rec.period = [$('#r-p0').value, $('#r-p1').value];
		rec.season = $('#r-season-on').checked ? [$('#r-s0').value.trim(), $('#r-s1').value.trim()] : null;
		rec.quality = {};
		for (const [k] of QUALITY) rec.quality[k] = qv($('#r-q-' + k).value);
		rec.quality.cloud_score = $('#r-cs-on').checked ? num('#r-cs') : null;
		const sc = num('#r-scene');
		if (sc !== null) rec.quality.scene_clouds = sc;
		rec.exclude_landcover = $$('#r-wc .chip.on').map(c => Number(c.dataset.v));
		rec.variables = S.recipe.variables;
		if (!rec.folder) throw new Error('Write the destination folder in inVars.');
		if (!/^[A-Za-z0-9_-]+$/.test(rec.folder)) throw new Error('The folder can only have letters (without accents), digits, - and _.');
		if (!rec.study.resolution) throw new Error('Write the pixel size of the study area.');
		if (rec.crs !== 'auto' && !rec.crs) throw new Error('Write the code of the coordinate system (e.g. EPSG:25830).');
		if (!rec.period[0] || !rec.period[1]) throw new Error('Write the period.');
		if (rec.season && !rec.season.every(s => /^\d\d-\d\d$/.test(s))) throw new Error('Write the season as MM-DD (e.g. 06-01).');
		if (!rec.variables.length) throw new Error('Add at least one variable.');
		return rec;
	}

	function syncEco() {
		const on = $('#r-eco-on').checked;
		for (const id of ['#r-eco-res', '#r-eco-crs', '#r-eco-drive']) $(id).disabled = !on || (id === '#r-eco-drive' && account() !== 'own');
	}
	$('#r-eco-on').addEventListener('change', syncEco);
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

	// World cover chips
	for (const [v, label] of [[10, 'Tree cover'], [20, 'Shrubland'], [30, 'Grassland'], [40, 'Cropland'], [50, 'Built-up'],
		[60, 'Bare / sparse vegetation'], [70, 'Snow and ice'], [80, 'Water'], [90, 'Herbaceous wetland'], [95, 'Mangroves'], [100, 'Moss and lichen']]) {
		const c = el('span', {class: 'chip', 'data-v': v}, `${v} ${label}`);
		c.addEventListener('click', () => c.classList.toggle('on'));
		$('#r-wc').append(c);
	}

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
					el('div', {class: 'var-src'}, k === 'catalog' ? 'catalog of MonitorEO' : k === 'image' ? 'Earth Engine image' : 'Earth Engine collection')),
				el('td', {class: 'small'}, src),
				el('td', {}, k === 'image' ? el('span', {class: 'faint'}, 'fixed image') : ms.map(m => el('span', {class: 'badge brand', style: {marginRight: '3px'}}, m))),
				el('td', {class: 'small'}, k === 'image' ? '—' : (TEMPORAL[v.temporal || 'period'] + (v.temporal && v.temporal !== 'period' ? ' · ' + (v.interannual || 'mean') : ''))),
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
				el('div', {class: 'small faint'}, `Catalog of MonitorEO (${(cat.source || {}).imported || ''}).`)),
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
				for (const k of ['period', 'season', 'quality', 'exclude_landcover', 'resample', 'aggregate']) if (existing && existing[k] !== undefined) v[k] = existing[k];
				return v;
			}}]});
	}

	async function eeDialog(existing) {
		const v0 = existing || {};
		const kindSel = el('select', {}, el('option', {value: 'image'}, 'Fixed image (e.g. a DEM)'), el('option', {value: 'collection'}, 'Image collection'));
		kindSel.value = existing ? kind(existing) : 'image';
		const id = el('input', {type: 'text', value: v0.image || v0.collection || '', placeholder: 'COPERNICUS/DEM/GLO30', style: {width: '100%'}});
		const band = el('input', {type: 'text', value: v0.band || '', placeholder: 'DEM'});
		const scale = el('input', {type: 'number', value: v0.scale ?? '', step: 'any', placeholder: '1'});
		const offset = el('input', {type: 'number', value: v0.offset ?? '', step: 'any', placeholder: '0'});
		const terrain = el('select', {}, el('option', {value: ''}, 'The value of the band'), el('option', {value: 'slope'}, 'Slope (degrees)'), el('option', {value: 'aspect'}, 'Aspect (degrees)'));
		terrain.value = v0.terrain || '';
		const mbox = metricChips(existing ? metrics(existing) : ['mean']);
		const tctl = temporalCtl(v0);
		const name = el('input', {type: 'text', value: v0.name || '', placeholder: 'automatic'});
		const rowsImg = [el('label', {}, 'Computation'), terrain];
		const rowsCol = [el('label', {}, 'Metrics'), mbox, el('label', {}, 'Temporal aggregation'), tctl.t, el('label', {}, 'Between years'), tctl.inter];
		const grid = el('div', {class: 'form-grid'}, el('label', {}, 'Type'), kindSel, el('label', {}, 'Earth Engine ID'), id,
			el('label', {}, 'Band'), band, el('label', {}, 'Scale and offset'), el('div', {class: 'row'}, scale, offset),
			...rowsImg, ...rowsCol, el('label', {}, 'File name'), name);
		const sync = () => {
			rowsImg.forEach(e => e.classList.toggle('hidden', kindSel.value !== 'image'));
			rowsCol.forEach(e => e.classList.toggle('hidden', kindSel.value !== 'collection'));
		};
		kindSel.addEventListener('change', sync); sync();
		return modal({title: existing ? 'Edit the variable' : 'Other collection or image of Earth Engine', body: el('div', {class: 'stack'},
			el('div', {class: 'small muted'}, 'Any dataset of the Earth Engine catalog (developers.google.com/earth-engine/datasets): '
				+ 'e.g. WORLDCLIM/V1/BIO band bio01 with scale 0.1 (mean annual temperature), or IDAHO_EPSCOR/TERRACLIMATE band pr.'), grid),
			buttons: [{label: 'Cancel', result: null}, {label: existing ? 'Save' : 'Add', cls: 'primary', value: () => {
				if (!id.value.trim() || !band.value.trim()) { toast('Write the ID and the band.', 'warn'); return undefined; }
				const v = {};
				if (name.value.trim()) v.name = name.value.trim();
				v[kindSel.value] = id.value.trim();
				v.band = band.value.trim();
				if (scale.value !== '') v.scale = Number(scale.value);
				if (offset.value !== '') v.offset = Number(offset.value);
				if (kindSel.value === 'image') { if (terrain.value) v.terrain = terrain.value; }
				else {
					const ms = chosenMetrics(mbox);
					if (!ms.length) { toast('Choose at least one metric.', 'warn'); return undefined; }
					v.metric = ms.length === 1 ? ms[0] : ms;
					if (tctl.t.value !== 'period') { v.temporal = tctl.t.value; v.interannual = tctl.inter.value; }
				}
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
		const opts = {dry_run: dry, only: $('#g-only').value || null, overwrite: $('#g-over').checked, use: $('#g-use').checked};
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
		},
		attach(job) { view.kind = job.kind === 'gee' ? 'download' : job.kind; view.start(job.id, job.title); },
		onShow() { sizeNote(); },
	};
})();

// 3. HFT mapping =====================================================================================================
pages.model = (() => {
	const S = {params: [], saved: {}, status: null};
	const view = new JobView($('#m-log'), $('#m-bar'), {
		onStart: () => { $('#m-stop').classList.remove('hidden'); $('#m-run').disabled = true; $('#dot-model').classList.remove('hidden'); $('#m-after').classList.add('hidden'); },
		onEnd: s => {
			$('#m-stop').classList.add('hidden'); $('#m-run').disabled = false; $('#dot-model').classList.add('hidden');
			refresh();
			if (s.state === 'done') {
				toast('Model finished.', 'ok');
				$('#m-after').classList.remove('hidden');
				pages.results.reload(s.result);
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
			f('Study area', st.selected_study_area || '—', st.study_area_file ? `${st.study_area_file} · IDs: ${st.ids ? st.ids.join(', ') : 'all'}` : 'choose it in Input variables'),
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

	async function refresh() {
		try {
			const st = await api('status');
			S.status = st;
			facts(st);
			steps(st);
			const sc = $('#p-SCENARIO');
			if (sc) sc.placeholder = st.scenario;
		} catch (e) { fail(e); }
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
		opt($('#x-id'), ids.map(i => [i, 'ID ' + i]), $('#x-id').value);
		$('#x-id').classList.toggle('hidden', ids.length < 2);
		$('#x-id-lab').classList.toggle('hidden', ids.length < 2);
		$('#x-meta').textContent = S.data.metadata || 'No processing_metadata.txt';
		await onId();
	}

	async function onId() {
		const rec = S.data.per_id[$('#x-id').value];
		if (!rec) return;
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
		$('#x-map-title').textContent = `HFT map${c.cut === 'results' ? '' : ' · ' + c.cut}`;
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
		S.outline = L.geoJSON(gj, {style: {color: '#111', weight: 2.2, fill: false}, interactive: false});
		if ($('#x-outline').checked) S.outline.addTo(S.map);
	}

	function fit() {
		S.map.invalidateSize();
		if (S.hft) S.map.fitBounds(S.hft.bounds, {padding: [16, 16]});
	}

	function tables(rec) {
		const box = $('#x-tables');
		box.textContent = '';
		for (const t of rec.tables) {
			if (t.file.startsWith('r2_') && t.file.endsWith('_hfts.csv')) continue; // long: shown in the R² plot
			const head = t.header;
			const tbl = el('table', {class: 't'});
			tbl.append(el('thead', {}, el('tr', {}, ...head.map(h => el('th', {title: h}, shortName(h))))));
			const sugCol = head.indexOf('suggested');
			const selCol = head.indexOf('selected');
			const body = el('tbody');
			for (const r of t.rows) {
				const hi = (sugCol >= 0 && r[sugCol] === '1') || (selCol >= 0 && r[selCol] === '1');
				body.append(el('tr', {class: hi ? 'sel' : ''}, ...r.map(v => el('td', {class: /^-?[\d.]+(e-?\d+)?$/.test(v) ? 'num' : ''}, /^-?\d+\.\d{4,}/.test(v) ? fmt(v, 3) : v))));
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

// 5. Similarity and landscape (not yet functional) ===================================================================
pages.similarity = {
	async onShow() {
		const list = await api('results_list').catch(() => []);
		const s = $('#s-sc');
		s.textContent = '';
		for (const x of list) for (const sc of x.scenarios) s.append(el('option', {}, `${x.study_area} / ${sc}`));
		if (!s.options.length) s.append(el('option', {}, '— no results —'));
		const sets = await api('sets').catch(() => ({sets: []}));
		const e = $('#s-set');
		e.textContent = '';
		for (const x of sets.sets.filter(x => x.ecoregion)) e.append(el('option', {}, `${x.label}/ecoregion · ${x.ecoregion} variables`));
		if (!e.options.length) e.append(el('option', {}, '— no ecoregion variables —'));
	},
};

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
		}
	} catch (e) { /* nothing */ }
	const p = store('page');
	if (p && pages[p] && p !== 'inputs') show(p);
})();
