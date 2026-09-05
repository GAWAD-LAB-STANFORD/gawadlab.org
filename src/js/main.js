// Gawad Lab site — small, dependency-free behaviours.
(function () {
  // Mobile navigation
  var toggle = document.querySelector('.nav-toggle');
  var links = document.querySelector('.nav-links');
  if (toggle && links) {
    toggle.addEventListener('click', function () {
      var open = links.classList.toggle('open');
      toggle.setAttribute('aria-expanded', open ? 'true' : 'false');
    });
    document.addEventListener('keydown', function (e) {
      if (e.key === 'Escape' && links.classList.contains('open')) {
        links.classList.remove('open'); toggle.setAttribute('aria-expanded', 'false');
      }
    });
  }

  // "Show all authors" toggles
  document.querySelectorAll('.pub-authors .link-btn').forEach(function (btn) {
    btn.addEventListener('click', function () {
      var p = btn.closest('.pub-authors');
      var expanded = p.classList.toggle('expanded');
      btn.textContent = expanded ? 'show fewer' : 'show all ' + btn.dataset.count + ' authors';
    });
  });

  // Lightbox for figures (a.zoom)
  var zooms = document.querySelectorAll('a.zoom');
  if (zooms.length) {
    var lb = document.createElement('div'); lb.className = 'lightbox'; lb.hidden = true;
    lb.innerHTML = '<button class="close" aria-label="Close">&times;</button><img alt=""><div class="cap"></div>';
    document.body.appendChild(lb);
    var lbImg = lb.querySelector('img'), lbCap = lb.querySelector('.cap');
    function closeLb() { lb.hidden = true; lbImg.src = ''; document.body.style.overflow = ''; }
    zooms.forEach(function (a) {
      a.addEventListener('click', function (e) {
        e.preventDefault();
        lbImg.src = a.getAttribute('href'); lbImg.alt = a.querySelector('img').alt;
        lbCap.textContent = a.querySelector('img').alt;
        lb.hidden = false; document.body.style.overflow = 'hidden';
      });
    });
    lb.addEventListener('click', closeLb);
    document.addEventListener('keydown', function (e) { if (e.key === 'Escape' && !lb.hidden) closeLb(); });
  }

  // Publications search + filter
  var search = document.getElementById('pub-search');
  if (!search) return;
  var pubs = Array.prototype.slice.call(document.querySelectorAll('.pub[data-search]'));
  var years = Array.prototype.slice.call(document.querySelectorAll('.pub-year[data-year]'));
  var chips = Array.prototype.slice.call(document.querySelectorAll('.chip[data-filter]'));
  var count = document.getElementById('pub-count');
  var empty = document.getElementById('pub-empty');
  var active = 'all';

  function apply() {
    var q = search.value.trim().toLowerCase();
    var shown = 0;
    pubs.forEach(function (li) {
      var okText = !q || li.dataset.search.indexOf(q) !== -1;
      var okFilter = active === 'all' || (active === 'highlight' && li.classList.contains('highlight')) ||
                     (active === 'preprint' && li.dataset.preprint === '1') ||
                     (active === 'lab' && li.dataset.lab === '1');
      var ok = okText && okFilter;
      li.hidden = !ok; if (ok) shown++;
    });
    years.forEach(function (h) {
      var any = pubs.some(function (li) { return !li.hidden && li.dataset.year === h.dataset.year; });
      h.hidden = !any;
    });
    if (count) count.textContent = shown + ' of ' + pubs.length + ' publications';
    if (empty) empty.hidden = shown !== 0;
  }
  search.addEventListener('input', apply);
  chips.forEach(function (c) {
    c.addEventListener('click', function () {
      active = c.dataset.filter;
      chips.forEach(function (o) { o.setAttribute('aria-pressed', o === c ? 'true' : 'false'); });
      apply();
    });
  });
  apply();
})();
