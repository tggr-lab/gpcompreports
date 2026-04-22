// landing_v2.js — landing-only behaviors: stat count-up + search teaser.
// Theme toggle, tooltips, collapsibles, reveals live in site_v2.js.

(function() {
  'use strict';

  var reduceMotion = matchMedia('(prefers-reduced-motion: reduce)').matches;

  // ------------------ Stat count-up on scroll-into-view ------------------
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

  // ------------------ Search teaser ------------------
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
      setTimeout(function() { resultsList.hidden = true; }, 120);
    });
    input.addEventListener('focus', function() {
      if (currentMatches.length) resultsList.hidden = false;
    });

    document.addEventListener('keydown', function(e) {
      if (e.key !== '/') return;
      var tag = (document.activeElement && document.activeElement.tagName) || '';
      if (tag === 'INPUT' || tag === 'TEXTAREA' || document.activeElement.isContentEditable) return;
      e.preventDefault();
      input.focus();
    });
  }
})();
