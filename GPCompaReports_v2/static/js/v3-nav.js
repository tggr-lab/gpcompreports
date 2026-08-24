/* v3-nav.js — sticky section nav for report pages.
   Builds itself from [data-section-title] so adding a section needs no JS change.
   Navigation is wayfinding, not content: this is always on. */
(function () {
  'use strict';
  var V3Nav = (window.V3Nav = window.V3Nav || {});

  function build() {
    var sections = document.querySelectorAll('.report-section[data-section-title]');
    if (sections.length < 2) return;

    var nav = document.createElement('nav');
    nav.className = 'v3-secnav';
    nav.setAttribute('aria-label', 'Report sections');
    var list = document.createElement('div');
    list.className = 'v3-secnav-list';
    nav.appendChild(list);

    var links = [];
    Array.prototype.forEach.call(sections, function (sec) {
      var a = document.createElement('a');
      a.href = '#' + sec.id;
      a.textContent = sec.dataset.sectionTitle;
      a.dataset.target = sec.id;
      list.appendChild(a);
      links.push(a);
    });

    var host = document.querySelector('.report-page .section-inner');
    var anchor = document.querySelector('.report-stats');
    if (!host || !anchor) return;
    anchor.parentNode.insertBefore(nav, anchor.nextSibling);

    if (!('IntersectionObserver' in window)) return;
    var seen = {};
    var obs = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) { seen[e.target.id] = e.isIntersecting; });
      var activeId = null;
      Array.prototype.forEach.call(sections, function (sec) {
        if (!activeId && seen[sec.id]) activeId = sec.id;
      });
      links.forEach(function (a) {
        a.classList.toggle('on', a.dataset.target === activeId);
      });
    }, { rootMargin: '-80px 0px -70% 0px' });
    Array.prototype.forEach.call(sections, function (s) { obs.observe(s); });
  }

  V3Nav.build = build;
  document.addEventListener('DOMContentLoaded', build);
})();
