// landing_v2.js — theme toggle, stat count-up, search teaser, scroll reveals.
// No frameworks. Runs on page load after initial theme boot (in <head>).

(function() {
  'use strict';

  var root = document.documentElement;
  var reduceMotion = matchMedia('(prefers-reduced-motion: reduce)').matches;

  // ------------------------------------------------------------------
  // Theme toggle — cycles auto -> light -> dark -> auto
  // ------------------------------------------------------------------
  var toggle = document.getElementById('theme-toggle');
  if (toggle) {
    function readStored() {
      try { return localStorage.getItem('gpcrsynth-theme'); }
      catch (e) { return null; }
    }
    function writeStored(v) {
      try {
        if (v) localStorage.setItem('gpcrsynth-theme', v);
        else   localStorage.removeItem('gpcrsynth-theme');
      } catch (e) {}
    }
    function applyMode(storedValue) {
      if (storedValue === 'light' || storedValue === 'dark') {
        root.setAttribute('data-color-mode', storedValue);
        root.removeAttribute('data-theme-auto');
      } else {
        var sysDark = matchMedia('(prefers-color-scheme: dark)').matches;
        root.setAttribute('data-color-mode', sysDark ? 'dark' : 'light');
        root.setAttribute('data-theme-auto', '');
      }
      // aria-pressed reflects "is theme overridden from auto"
      toggle.setAttribute('aria-pressed', storedValue ? 'true' : 'false');
      toggle.title = storedValue
        ? 'Theme: ' + storedValue + ' (click to cycle)'
        : 'Theme: auto (click to cycle)';
    }

    // Sync to whatever boot script established
    applyMode(readStored());

    toggle.addEventListener('click', function() {
      var current = readStored();  // null | 'light' | 'dark'
      var next;
      if (current === null)       next = 'light';
      else if (current === 'light') next = 'dark';
      else                         next = null;  // back to auto
      writeStored(next);
      applyMode(next);
    });

    // Follow system changes while in auto mode
    try {
      matchMedia('(prefers-color-scheme: dark)').addEventListener('change', function() {
        if (readStored() === null) applyMode(null);
      });
    } catch (e) { /* older Safari */ }
  }

  // ------------------------------------------------------------------
  // Stat count-up on scroll-into-view
  // ------------------------------------------------------------------
  var statNums = document.querySelectorAll('.stat-num[data-target]');
  function easeOutQuart(t) { return 1 - Math.pow(1 - t, 4); }
  function animateCount(el) {
    var target = parseFloat(el.getAttribute('data-target'));
    var suffix = el.getAttribute('data-suffix') || '';
    if (!isFinite(target)) return;
    if (reduceMotion) { el.textContent = target + suffix; return; }
    var duration = 1200;
    var start = performance.now();
    function tick(now) {
      var t = Math.min(1, (now - start) / duration);
      var value = Math.round(target * easeOutQuart(t));
      el.textContent = value + suffix;
      if (t < 1) requestAnimationFrame(tick);
    }
    requestAnimationFrame(tick);
  }
  if ('IntersectionObserver' in window && statNums.length) {
    var statObs = new IntersectionObserver(function(entries) {
      entries.forEach(function(entry) {
        if (entry.isIntersecting) {
          animateCount(entry.target);
          statObs.unobserve(entry.target);
        }
      });
    }, { threshold: 0.5 });
    statNums.forEach(function(el) { statObs.observe(el); });
  } else {
    statNums.forEach(function(el) {
      el.textContent = el.getAttribute('data-target') + (el.getAttribute('data-suffix') || '');
    });
  }

  // ------------------------------------------------------------------
  // Scroll reveals
  // ------------------------------------------------------------------
  var reveals = document.querySelectorAll('.reveal');
  if ('IntersectionObserver' in window && reveals.length && !reduceMotion) {
    var revealObs = new IntersectionObserver(function(entries) {
      entries.forEach(function(entry) {
        if (entry.isIntersecting) {
          entry.target.classList.add('is-visible');
          revealObs.unobserve(entry.target);
        }
      });
    }, { threshold: 0.08, rootMargin: '0px 0px -40px 0px' });
    reveals.forEach(function(el) { revealObs.observe(el); });
  } else {
    reveals.forEach(function(el) { el.classList.add('is-visible'); });
  }

  // ------------------------------------------------------------------
  // Search teaser
  // ------------------------------------------------------------------
  var input = document.getElementById('search-input');
  var resultsList = document.getElementById('search-results');
  var dataScript = document.getElementById('gpcr-search-data');
  var records = [];
  try { records = dataScript ? JSON.parse(dataScript.textContent) : []; }
  catch (e) { records = []; }

  if (input && resultsList && records.length) {
    var selected = -1;
    var currentMatches = [];

    function rankMatches(query) {
      if (!query) return [];
      var q = query.toLowerCase();
      var scored = [];
      for (var i = 0; i < records.length; i++) {
        var r = records[i];
        var gene = (r.gene_name || '').toLowerCase();
        var uniprot = (r.uniprot_name || '').toLowerCase();
        var id = (r.gpcr_id || '').toLowerCase();
        var fam = (r.receptor_family || '').toLowerCase();
        var score = -1;
        if (gene === q || uniprot === q) score = 0;
        else if (gene.startsWith(q) || uniprot.startsWith(q)) score = 1;
        else if (gene.indexOf(q) !== -1 || uniprot.indexOf(q) !== -1) score = 2;
        else if (id.indexOf(q) !== -1) score = 3;
        else if (fam.indexOf(q) !== -1) score = 4;
        if (score >= 0) scored.push({ score: score, record: r });
      }
      scored.sort(function(a, b) { return a.score - b.score; });
      return scored.slice(0, 6).map(function(s) { return s.record; });
    }

    function renderResults() {
      resultsList.innerHTML = '';
      if (!currentMatches.length) { resultsList.hidden = true; return; }
      currentMatches.forEach(function(r, idx) {
        var li = document.createElement('li');
        li.setAttribute('role', 'option');
        li.dataset.id = r.gpcr_id;
        li.setAttribute('aria-selected', idx === selected ? 'true' : 'false');
        li.innerHTML =
          '<span class="result-gene">' + escapeHtml(r.gene_name || r.uniprot_name) + '</span>' +
          '<span class="result-name">' + escapeHtml(r.uniprot_name) + '</span>' +
          '<span class="result-fam">' + escapeHtml(r.receptor_family || '') + '</span>';
        li.addEventListener('mouseenter', function() { selected = idx; updateSelection(); });
        li.addEventListener('mousedown', function(e) { e.preventDefault(); navigateTo(r.gpcr_id); });
        resultsList.appendChild(li);
      });
      resultsList.hidden = false;
    }

    function updateSelection() {
      Array.prototype.forEach.call(resultsList.children, function(li, idx) {
        li.setAttribute('aria-selected', idx === selected ? 'true' : 'false');
      });
    }

    function navigateTo(gpcrId) {
      if (gpcrId) window.location.href = 'reports/' + gpcrId + '.html';
    }

    function escapeHtml(s) {
      return String(s).replace(/[&<>"']/g, function(c) {
        return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c];
      });
    }

    input.addEventListener('input', function() {
      currentMatches = rankMatches(input.value.trim());
      selected = currentMatches.length ? 0 : -1;
      renderResults();
    });

    input.addEventListener('keydown', function(e) {
      if (!currentMatches.length && e.key !== 'Escape') return;
      if (e.key === 'ArrowDown') {
        e.preventDefault();
        selected = (selected + 1) % currentMatches.length;
        updateSelection();
      } else if (e.key === 'ArrowUp') {
        e.preventDefault();
        selected = (selected - 1 + currentMatches.length) % currentMatches.length;
        updateSelection();
      } else if (e.key === 'Enter') {
        e.preventDefault();
        if (selected >= 0 && currentMatches[selected]) navigateTo(currentMatches[selected].gpcr_id);
      } else if (e.key === 'Escape') {
        currentMatches = [];
        selected = -1;
        renderResults();
        input.blur();
      }
    });

    input.addEventListener('blur', function() {
      // Small delay so mousedown on a result can fire first
      setTimeout(function() { resultsList.hidden = true; }, 120);
    });
    input.addEventListener('focus', function() {
      if (currentMatches.length) resultsList.hidden = false;
    });

    // "/" focuses search (unless user is already typing in a field)
    document.addEventListener('keydown', function(e) {
      if (e.key !== '/') return;
      var tag = (document.activeElement && document.activeElement.tagName) || '';
      if (tag === 'INPUT' || tag === 'TEXTAREA' || document.activeElement.isContentEditable) return;
      e.preventDefault();
      input.focus();
    });
  }
})();
