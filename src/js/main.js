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

  // Support modal (replaces the bare mailto for the header/footer buttons)
  var modal = document.getElementById('support-modal');
  if (modal) {
    var openers = document.querySelectorAll('a.support-link');
    var lastFocus = null;
    function openModal(e) { e.preventDefault(); lastFocus = document.activeElement; modal.hidden = false; document.body.style.overflow = 'hidden'; modal.querySelector('.modal-close').focus(); }
    function closeModal() { modal.hidden = true; document.body.style.overflow = ''; if (lastFocus) lastFocus.focus(); }
    openers.forEach(function (a) { a.addEventListener('click', openModal); });
    modal.querySelector('.modal-close').addEventListener('click', closeModal);
    modal.addEventListener('click', function (e) { if (e.target === modal) closeModal(); });
    document.addEventListener('keydown', function (e) { if (e.key === 'Escape' && !modal.hidden) closeModal(); });
    var copyBtn = document.getElementById('copy-email');
    if (copyBtn) copyBtn.addEventListener('click', function () {
      var addr = document.getElementById('support-address').textContent.trim();
      var done = function () { copyBtn.textContent = 'Copied'; setTimeout(function () { copyBtn.textContent = 'Copy email'; }, 2000); };
      if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(addr).then(done, done);
      else { var r = document.createRange(); r.selectNodeContents(document.getElementById('support-address')); var sel = window.getSelection(); sel.removeAllRanges(); sel.addRange(r); try { document.execCommand('copy'); } catch (err) {} done(); }
    });
  }

  // Header compacts after scrolling
  var header = document.querySelector('.site-header');
  if (header) {
    var onScroll = function () { header.classList.toggle('scrolled', window.scrollY > 60); };
    window.addEventListener('scroll', onScroll, { passive: true }); onScroll();
  }

  // Scroll reveal
  var reduce = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  if (!reduce && 'IntersectionObserver' in window) {
    var singles = document.querySelectorAll('.section-head, .feature, .figure-wide, .prose, .contact-card, .cta-band .container, .growing .container, .photo-band .container, .photo-duo, .pub-year, .page-hero .container');
    var groups = document.querySelectorAll('.pillars, .research-index, .pub-compact, .timeline, .people-grid, .alumni-grid, .info-grid, .why-grid, .resource-groups, .openings, .mascots, .pub-list, .translate-grid');
    singles.forEach(function (el) { el.classList.add('reveal'); });
    groups.forEach(function (el) { el.classList.add('reveal-stagger'); });
    var io = new IntersectionObserver(function (entries) {
      entries.forEach(function (en) { if (en.isIntersecting) { en.target.classList.add('in'); io.unobserve(en.target); } });
    }, { rootMargin: '0px 0px -8% 0px', threshold: 0.08 });
    document.querySelectorAll('.reveal, .reveal-stagger').forEach(function (el) { io.observe(el); });
    // Fallbacks so content can never stay hidden: check positions on load and on scroll,
    // and reveal everything if the observer has not fired within a few seconds.
    var revealByPosition = function () {
      document.querySelectorAll('.reveal:not(.in), .reveal-stagger:not(.in)').forEach(function (el) {
        var r = el.getBoundingClientRect(); if (r.top < window.innerHeight * 0.95) el.classList.add('in');
      });
    };
    setTimeout(revealByPosition, 50);
    var ticking = false;
    window.addEventListener('scroll', function () { if (!ticking) { ticking = true; requestAnimationFrame(function () { revealByPosition(); ticking = false; }); } }, { passive: true });
    setTimeout(function () { if (document.visibilityState !== 'visible') document.querySelectorAll('.reveal, .reveal-stagger').forEach(function (el) { el.classList.add('in'); }); }, 4000);
  }

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
