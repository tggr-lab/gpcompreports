/* v3-deeplink.js — put snake-plot state in the URL hash so a link can point at a view.
   Grammar: #view=cfr&min=0.5&max=3&dir=active&sec=sec-changes
   Reading the hash changes what is displayed but adds nothing to the page. */
(function () {
  'use strict';

  function parseHash() {
    var raw = (window.location.hash || '').replace(/^#/, '');
    var out = {};
    if (!raw) return out;
    raw.split('&').forEach(function (pair) {
      var kv = pair.split('=');
      if (kv.length === 2 && kv[0]) out[decodeURIComponent(kv[0])] = decodeURIComponent(kv[1]);
    });
    return out;
  }

  function writeHash(state) {
    var parts = [];
    if (state.view && state.view !== 'delta_rrcs') parts.push('view=' + state.view);
    var b = window.SnakeAPI.bounds();
    if (typeof state.min === 'number' && state.min > b.lo) parts.push('min=' + state.min.toFixed(2));
    if (typeof state.max === 'number' && state.max < b.hi) parts.push('max=' + state.max.toFixed(2));
    if (state.dir && state.dir !== 'both') parts.push('dir=' + state.dir);
    var hash = parts.length ? '#' + parts.join('&') : '';
    if (hash !== window.location.hash) {
      history.replaceState(null, '', window.location.pathname + window.location.search + hash);
    }
  }

  // `sec` just scrolls to a report section and has no dependency on the snake
  // plot, so it is applied on its own and must work even when SnakeAPI never
  // loads (e.g. a page whose snake plot build failed).
  function applySectionHash() {
    var h = parseHash();
    if (h.sec) {
      var target = document.getElementById(h.sec);
      if (target) target.scrollIntoView();
    }
  }

  function applySnakeHash() {
    var h = parseHash();
    if (!window.SnakeAPI) return;
    if (h.view) window.SnakeAPI.setView(h.view);
    if (h.min || h.max) {
      window.SnakeAPI.setLinksVisible(true);
      var cb = document.getElementById('snake-links-cb');
      if (cb) cb.checked = true;
      window.SnakeAPI.setMagRange(
        h.min ? parseFloat(h.min) : undefined,
        h.max ? parseFloat(h.max) : undefined);
    }
    if (h.dir) window.SnakeAPI.setDirection(h.dir);
  }

  function applyHash() {
    applySnakeHash();
    applySectionHash();
  }

  function watch() {
    var sync = function () { writeHash(window.SnakeAPI.getState()); };
    document.querySelectorAll('.snake-btn').forEach(function (b) {
      b.addEventListener('click', function () { setTimeout(sync, 0); });
    });
    document.querySelectorAll('.snake-dir-btn').forEach(function (b) {
      b.addEventListener('click', function () { setTimeout(sync, 0); });
    });
    ['snake-mag-min', 'snake-mag-max'].forEach(function (id) {
      var el = document.getElementById(id);
      if (el) el.addEventListener('change', sync);
    });
  }

  function startSnake() {
    applySnakeHash();
    watch();
  }

  // `sec` has no SnakeAPI dependency, so it applies immediately and
  // unconditionally, once, regardless of whether a snake plot exists on this
  // page.
  applySectionHash();

  // Ordering: the snake plot's code is an INLINE script, so it runs at parse
  // time. This module is deferred, so it runs after parsing, by which point
  // SnakeAPI already exists and `snake-ready` has already fired. Check for the
  // object first and only fall back to the event, otherwise a listener attached
  // here would wait for an event that already happened and never apply a hash.
  if (window.SnakeAPI) {
    startSnake();
  } else {
    document.addEventListener('snake-ready', startSnake);
  }

  // Registered exactly once, unconditionally, so both halves of the grammar
  // stay live on hashchange whether or not a snake plot exists on this page.
  window.addEventListener('hashchange', applyHash);
})();
