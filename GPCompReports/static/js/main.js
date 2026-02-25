/* GPCompReports — Client-side search, filter, sort, pagination + Material interactions */

(function() {
    'use strict';

    let gpcrData = [];
    let filteredData = [];
    let currentPage = 1;
    let perPage = 25;
    let sortCol = 'sum_abs_delta';
    let sortDir = 'desc';

    // Initialize when DOM is ready
    document.addEventListener('DOMContentLoaded', function() {
        const tableEl = document.getElementById('gpcr-table-body');
        if (!tableEl) return;

        loadData();
    });

    function loadData() {
        // Use inline data embedded in the page (works with file:// protocol)
        if (window.GPCR_DATA) {
            gpcrData = window.GPCR_DATA;
            filteredData = gpcrData.slice();
            populateFilters();
            sortData();
            render();
            bindEvents();
            return;
        }
        // Fallback to fetch
        fetch('../data/gpcr_index.json')
            .then(r => r.json())
            .then(data => {
                gpcrData = data;
                filteredData = data.slice();
                populateFilters();
                sortData();
                render();
                bindEvents();
            })
            .catch(e => console.error('Failed to load GPCR data:', e));
    }

    function populateFilters() {
        const ligandSelect = document.getElementById('filter-ligand');
        const familySelect = document.getElementById('filter-family');
        if (!ligandSelect || !familySelect) return;

        const ligands = [...new Set(gpcrData.map(d => d.ligand_type))].sort();
        const families = [...new Set(gpcrData.map(d => d.receptor_family))].sort();

        ligands.forEach(l => {
            const opt = document.createElement('option');
            opt.value = l; opt.textContent = l;
            ligandSelect.appendChild(opt);
        });

        families.forEach(f => {
            const opt = document.createElement('option');
            opt.value = f; opt.textContent = f;
            familySelect.appendChild(opt);
        });
    }

    function bindEvents() {
        const searchInput = document.getElementById('search-input');
        const ligandSelect = document.getElementById('filter-ligand');
        const familySelect = document.getElementById('filter-family');
        const perPageSelect = document.getElementById('per-page');

        if (searchInput) searchInput.addEventListener('input', applyFilters);
        if (ligandSelect) ligandSelect.addEventListener('change', applyFilters);
        if (familySelect) familySelect.addEventListener('change', applyFilters);
        if (perPageSelect) perPageSelect.addEventListener('change', function() {
            perPage = parseInt(this.value);
            currentPage = 1;
            render();
        });

        // Column sort headers
        document.querySelectorAll('[data-sort]').forEach(th => {
            th.addEventListener('click', function() {
                const col = this.dataset.sort;
                if (sortCol === col) {
                    sortDir = sortDir === 'asc' ? 'desc' : 'asc';
                } else {
                    sortCol = col;
                    sortDir = 'desc';
                }
                sortData();
                render();
                updateSortIndicators();
            });
        });
    }

    function applyFilters() {
        const query = (document.getElementById('search-input')?.value || '').toLowerCase();
        const ligand = document.getElementById('filter-ligand')?.value || '';
        const family = document.getElementById('filter-family')?.value || '';

        filteredData = gpcrData.filter(d => {
            const matchSearch = !query ||
                d.uniprot_name.toLowerCase().includes(query) ||
                d.gene_name.toLowerCase().includes(query) ||
                d.receptor_family.toLowerCase().includes(query) ||
                d.gpcr_id.toLowerCase().includes(query);
            const matchLigand = !ligand || d.ligand_type === ligand;
            const matchFamily = !family || d.receptor_family === family;
            return matchSearch && matchLigand && matchFamily;
        });

        currentPage = 1;
        sortData();
        render();
    }

    function sortData() {
        filteredData.sort((a, b) => {
            let va = a[sortCol], vb = b[sortCol];
            if (typeof va === 'string') {
                va = va.toLowerCase(); vb = vb.toLowerCase();
                return sortDir === 'asc' ? va.localeCompare(vb) : vb.localeCompare(va);
            }
            return sortDir === 'asc' ? va - vb : vb - va;
        });
    }

    function render() {
        const tbody = document.getElementById('gpcr-table-body');
        if (!tbody) return;

        const totalPages = Math.ceil(filteredData.length / perPage);
        const start = (currentPage - 1) * perPage;
        const pageData = filteredData.slice(start, start + perPage);

        tbody.innerHTML = pageData.map((d, i) => `
            <tr onclick="window.location='../reports/${d.gpcr_id}.html'" class="clickable-row">
                <td>${start + i + 1}</td>
                <td><a href="../reports/${d.gpcr_id}.html">${d.uniprot_name}</a></td>
                <td>${d.gene_name}</td>
                <td>${d.receptor_family}</td>
                <td><span class="badge badge-teal">${d.ligand_type}</span></td>
                <td>${d.total_contacts}</td>
                <td>${d.significant_changes}</td>
                <td><strong>${d.sum_abs_delta.toFixed(1)}</strong></td>
                <td>${d.variants_found}</td>
            </tr>
        `).join('');

        renderPagination(totalPages);
        updateResultCount();
    }

    function renderPagination(totalPages) {
        const container = document.getElementById('pagination');
        if (!container) return;

        let html = '';
        html += `<button onclick="GPComp.goPage(${currentPage - 1})" ${currentPage <= 1 ? 'disabled' : ''}>← Prev</button>`;

        const maxButtons = 7;
        let startPage = Math.max(1, currentPage - Math.floor(maxButtons / 2));
        let endPage = Math.min(totalPages, startPage + maxButtons - 1);
        if (endPage - startPage < maxButtons - 1) {
            startPage = Math.max(1, endPage - maxButtons + 1);
        }

        if (startPage > 1) {
            html += `<button onclick="GPComp.goPage(1)">1</button>`;
            if (startPage > 2) html += `<span class="page-info">...</span>`;
        }

        for (let p = startPage; p <= endPage; p++) {
            html += `<button onclick="GPComp.goPage(${p})" class="${p === currentPage ? 'active' : ''}">${p}</button>`;
        }

        if (endPage < totalPages) {
            if (endPage < totalPages - 1) html += `<span class="page-info">...</span>`;
            html += `<button onclick="GPComp.goPage(${totalPages})">${totalPages}</button>`;
        }

        html += `<button onclick="GPComp.goPage(${currentPage + 1})" ${currentPage >= totalPages ? 'disabled' : ''}>Next →</button>`;

        container.innerHTML = html;
    }

    function updateResultCount() {
        const el = document.getElementById('result-count');
        if (el) {
            el.textContent = `${filteredData.length} of ${gpcrData.length} GPCRs`;
        }
    }

    function updateSortIndicators() {
        document.querySelectorAll('[data-sort]').forEach(th => {
            const icon = th.querySelector('.sort-icon');
            if (!icon) return;
            if (th.dataset.sort === sortCol) {
                th.classList.add('sorted');
                icon.textContent = sortDir === 'asc' ? ' ▲' : ' ▼';
            } else {
                th.classList.remove('sorted');
                icon.textContent = ' ⇅';
            }
        });
    }

    // Expose for onclick handlers
    window.GPComp = {
        goPage: function(page) {
            const totalPages = Math.ceil(filteredData.length / perPage);
            if (page < 1 || page > totalPages) return;
            currentPage = page;
            render();
            document.getElementById('gpcr-table')?.scrollIntoView({ behavior: 'smooth', block: 'start' });
        }
    };
})();

/* ===== Report Page Table Interactivity ===== */

// Pagination state for report-page tables
const paginationState = {};

function initPagination(tableId, pageSize) {
    pageSize = pageSize || 25;
    if (!paginationState[tableId]) {
        paginationState[tableId] = { currentPage: 1, pageSize: pageSize, totalRows: 0 };
    }
    updatePagination(tableId);
}

function updatePagination(tableId) {
    var table = document.getElementById(tableId);
    if (!table) return;
    var tbody = table.getElementsByTagName('tbody')[0];
    if (!tbody) return;
    var rows = Array.from(tbody.getElementsByTagName('tr'));

    var visibleRows = rows.filter(function(r) { return r.style.display !== 'none'; });
    var state = paginationState[tableId];
    state.totalRows = visibleRows.length;
    var totalPages = Math.ceil(state.totalRows / state.pageSize);

    rows.forEach(function(r) { r.classList.add('pagination-hidden'); });

    var start = (state.currentPage - 1) * state.pageSize;
    var end = start + state.pageSize;
    visibleRows.slice(start, end).forEach(function(r) { r.classList.remove('pagination-hidden'); });

    var paginationDiv = document.getElementById('pagination-' + tableId);
    if (paginationDiv && totalPages > 1) {
        var pageInfo = 'Page ' + state.currentPage + ' of ' + totalPages + ' (' + state.totalRows + ' rows)';
        var prevDis = state.currentPage <= 1 ? 'disabled' : '';
        var nextDis = state.currentPage >= totalPages ? 'disabled' : '';
        paginationDiv.innerHTML =
            '<div class="pagination-controls">' +
            '<button onclick="changePage(\'' + tableId + '\',' + (state.currentPage - 1) + ')" ' + prevDis + ' class="btn">← Prev</button>' +
            '<span class="page-info">' + pageInfo + '</span>' +
            '<button onclick="changePage(\'' + tableId + '\',' + (state.currentPage + 1) + ')" ' + nextDis + ' class="btn">Next →</button>' +
            '<select onchange="changePageSize(\'' + tableId + '\',this.value)" class="page-size-select">' +
            '<option value="25"' + (state.pageSize === 25 ? ' selected' : '') + '>25 per page</option>' +
            '<option value="50"' + (state.pageSize === 50 ? ' selected' : '') + '>50 per page</option>' +
            '<option value="100"' + (state.pageSize === 100 ? ' selected' : '') + '>100 per page</option>' +
            '<option value="99999"' + (state.pageSize === 99999 ? ' selected' : '') + '>All</option>' +
            '</select></div>';
    } else if (paginationDiv) {
        paginationDiv.innerHTML = '<div class="pagination-controls"><span class="page-info">' + state.totalRows + ' rows</span></div>';
    }
}

function changePage(tableId, newPage) {
    var state = paginationState[tableId];
    var totalPages = Math.ceil(state.totalRows / state.pageSize);
    if (newPage >= 1 && newPage <= totalPages) {
        state.currentPage = newPage;
        updatePagination(tableId);
    }
}

function changePageSize(tableId, newSize) {
    var state = paginationState[tableId];
    state.pageSize = parseInt(newSize);
    state.currentPage = 1;
    updatePagination(tableId);
}

function searchTable(tableId) {
    var input = document.getElementById('search-' + tableId);
    if (!input) return;
    var filter = input.value.toUpperCase();
    var table = document.getElementById(tableId);
    if (!table) return;
    var tr = table.getElementsByTagName('tr');

    for (var i = 1; i < tr.length; i++) {
        var visible = false;
        var td = tr[i].getElementsByTagName('td');
        for (var j = 0; j < td.length; j++) {
            if (td[j]) {
                var text = td[j].textContent || td[j].innerText;
                if (text.toUpperCase().indexOf(filter) > -1) {
                    visible = true;
                    break;
                }
            }
        }
        tr[i].style.display = visible ? '' : 'none';
    }

    if (paginationState[tableId]) {
        paginationState[tableId].currentPage = 1;
        updatePagination(tableId);
    }
}

function sortTable(tableId, columnIndex) {
    var table = document.getElementById(tableId);
    if (!table) return;
    var th = table.getElementsByTagName('th')[columnIndex];
    var tbody = table.getElementsByTagName('tbody')[0];
    var rows = Array.from(tbody.getElementsByTagName('tr'));

    var isAsc = th.classList.contains('asc');
    Array.from(table.getElementsByTagName('th')).forEach(function(h) {
        h.classList.remove('asc', 'desc');
    });
    th.classList.add(isAsc ? 'desc' : 'asc');

    rows.sort(function(a, b) {
        var aVal = a.getElementsByTagName('td')[columnIndex].textContent.trim();
        var bVal = b.getElementsByTagName('td')[columnIndex].textContent.trim();
        var aNum = parseFloat(aVal.replace(/[^0-9.eE\-+]/g, ''));
        var bNum = parseFloat(bVal.replace(/[^0-9.eE\-+]/g, ''));
        if (!isNaN(aNum) && !isNaN(bNum)) {
            return isAsc ? bNum - aNum : aNum - bNum;
        }
        return isAsc ? bVal.localeCompare(aVal) : aVal.localeCompare(bVal);
    });

    rows.forEach(function(r) { tbody.appendChild(r); });

    if (paginationState[tableId]) {
        paginationState[tableId].currentPage = 1;
        updatePagination(tableId);
    }
}

function exportTableToCSV(tableId, filename) {
    var table = document.getElementById(tableId);
    if (!table) return;
    var rows = table.getElementsByTagName('tr');
    var csv = [];

    for (var i = 0; i < rows.length; i++) {
        if (rows[i].style.display === 'none') continue;
        var row = [];
        var cols = rows[i].querySelectorAll('td, th');
        for (var j = 0; j < cols.length; j++) {
            var text = cols[j].innerText.replace(/"/g, '""');
            row.push('"' + text + '"');
        }
        csv.push(row.join(','));
    }

    var blob = new Blob([csv.join('\n')], { type: 'text/csv' });
    var link = document.createElement('a');
    link.download = filename;
    link.href = window.URL.createObjectURL(blob);
    link.style.display = 'none';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
}

/* ===== Collapsible Cards ===== */

function toggleCollapse(cardId) {
    var body = document.getElementById(cardId + '-body');
    var icon = document.getElementById(cardId + '-icon');
    if (!body) return;
    body.classList.toggle('collapsed');
    if (icon) {
        icon.textContent = body.classList.contains('collapsed') ? '\u25B6' : '\u25BC';
    }
}

/* ===== Tooltip Smart Repositioning ===== */

function initTooltips() {
    var tooltipEls = document.querySelectorAll('[data-tooltip]');
    tooltipEls.forEach(function(el) {
        el.addEventListener('mouseenter', function() {
            var rect = el.getBoundingClientRect();
            var viewW = window.innerWidth;

            // Reset positioning classes
            el.classList.remove('tooltip-bottom', 'tooltip-left', 'tooltip-right');

            // Flip to bottom if too close to top
            if (rect.top < 80) {
                el.classList.add('tooltip-bottom');
            }

            // Shift right if too close to left edge
            if (rect.left < 160) {
                el.classList.add('tooltip-right');
            }

            // Shift left if too close to right edge
            if (viewW - rect.right < 160) {
                el.classList.add('tooltip-left');
            }
        });
    });
}

// Auto-init pagination for report tables and tooltips on DOMContentLoaded
document.addEventListener('DOMContentLoaded', function() {
    // Initialize Materialize components (sidenav, etc.)
    if (typeof M !== 'undefined') M.AutoInit();

    var tables = document.querySelectorAll('table[id]');
    tables.forEach(function(table) {
        var tableId = table.id;
        var pDiv = document.getElementById('pagination-' + tableId);
        if (pDiv) {
            initPagination(tableId, 25);
        }
    });

    // Initialize tooltip repositioning
    initTooltips();
});
