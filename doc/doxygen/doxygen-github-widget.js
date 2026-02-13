/**
 * SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
 * SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
 */

/**
 * Doxygen header GitHub widget
 *
 * What it does:
 *  - Adds a small "OWNER/REPO + stars + forks" widget in the top bar (near search).
 *  - Fetches live counts from the GitHub REST API.
 *  - Caches the counts in localStorage to reduce API calls / avoid rate limits.
 *  - Falls back gracefully if offline or blocked (shows cached values or placeholders).
 *
 * Placement:
 *  - Inserts the widget to the LEFT of the search box when #MSearchBox exists.
 *  - Otherwise, inserts into the main navigation menu (#main-menu).
 */
(function () {
    // ---------------------------------------------------------------------------
    // Configuration
    // ---------------------------------------------------------------------------
    const OWNER = "chaos-polymtl";
    const REPO  = "lethe";

    // GitHub links + API endpoint
    const REPO_URL = `https://github.com/${OWNER}/${REPO}`;
    const API_REPO = `https://api.github.com/repos/${OWNER}/${REPO}`;

    // localStorage cache settings
    const CACHE_KEY = `gh_widget:${OWNER}/${REPO}`;
    const CACHE_TTL_MS = 60 * 60 * 1000; // 1 hour

    // ---------------------------------------------------------------------------
    // Helpers
    // ---------------------------------------------------------------------------

    /**
     * Format a number with locale separators (e.g., 12345 -> "12,345").
     * If the value is missing, return a placeholder dash.
     */
    function formatCount(n) {
        if (typeof n !== "number") return "—";
        return n.toLocaleString(undefined);
    }

    /**
     * Read cached stats from localStorage, if present and still fresh.
     * Cache shape: { ts: <timestamp_ms>, data: { stargazers_count, forks_count } }
     */
    function loadCache() {
        try {
            const raw = localStorage.getItem(CACHE_KEY);
            if (!raw) return null;

            const obj = JSON.parse(raw);
            if (!obj?.ts || !obj?.data) return null;

            // Expire cache after TTL
            if (Date.now() - obj.ts > CACHE_TTL_MS) return null;

            return obj.data;
        } catch {
            // localStorage might be blocked or JSON might be corrupted
            return null;
        }
    }

    /**
     * Save stats to localStorage with a timestamp.
     */
    function saveCache(data) {
        try {
            localStorage.setItem(CACHE_KEY, JSON.stringify({ ts: Date.now(), data }));
        } catch {
            // localStorage might be unavailable (private mode, strict policies, etc.)
        }
    }

    /**
     * Find where to insert the widget.
     * Priority order:
     *  1. Before the search box (#MSearchBox) if it exists
     *  2. Inside the main menu (#main-menu) as a list item
     *  3. Inside #top header
     *  4. Fallback to body
     */
    function findInsertTarget() {
        // Try to find search box first (preferred placement)
        const search = document.getElementById("MSearchBox");
        if (search?.parentNode) {
            return { container: search.parentNode, before: search, type: 'search' };
        }

        // Try to find main menu (common on all Doxygen pages)
        const mainMenu = document.getElementById("main-menu");
        if (mainMenu) {
            return { container: mainMenu, before: null, type: 'menu' };
        }

        // Fallback to top header
        const top = document.getElementById("top");
        if (top) {
            return { container: top, before: null, type: 'header' };
        }

        // Last resort: body
        return { container: document.body, before: null, type: 'body' };
    }

    // ---------------------------------------------------------------------------
    // Icons (SVG using `currentColor` so CSS can control the color)
    // ---------------------------------------------------------------------------
    const ICONS = {
        github: `
      <svg viewBox="0 0 16 16" aria-hidden="true">
        <path fill="currentColor" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38
        0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52
        -.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78
        -.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82
        .64-.18 1.32-.27 2-.27s1.36.09 2 .27c1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12
        .51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93
        -.01 2.2 0 .21.15.46.55.38A8.01 8.01 0 0 0 16 8c0-4.42-3.58-8-8-8Z"></path>
      </svg>`,
        star: `
      <svg viewBox="0 0 24 24" aria-hidden="true">
        <path fill="currentColor" d="M12 17.27 18.18 21l-1.64-7.03L22 9.24l-7.19-.61L12 2 9.19 8.63 2 9.24l5.46 4.73L5.82 21z"/>
      </svg>`,
        fork: `
      <svg viewBox="0 0 24 24" aria-hidden="true">
        <path fill="currentColor" d="M7 4a3 3 0 1 0 2 2.83V8c0 1.66 1.34 3 3 3h1v6.17A3 3 0 1 0 15 20.17V11h1c1.66 0 3-1.34 3-3V6.83A3 3 0 1 0 17 4a3 3 0 0 0 0 5.83V8c0 .55-.45 1-1 1h-4c-.55 0-1-.45-1-1V9.83A3 3 0 0 0 7 4Z"/>
      </svg>`
    };

    // ---------------------------------------------------------------------------
    // Document Object Model (DOM) building / insertion
    // ---------------------------------------------------------------------------

    /**
     * Build the widget DOM node from the latest data.
     * `data` shape: { stargazers_count?: number, forks_count?: number }
     * `isMenuItem` determines if we wrap in <li> for menu insertion
     */
    function buildWidget(data, isMenuItem = false) {
        const repoName = `${OWNER}/${REPO}`;
        const stars    = formatCount(data?.stargazers_count);
        const forks    = formatCount(data?.forks_count);

        // Build the widget content
        const widgetHTML = `
      <div class="lgw-row1">
        <span class="lgw-icon lgw-gh">${ICONS.github}</span>
        <span class="lgw-name">${repoName}</span>
      </div>

      <div class="lgw-row2">
        <span class="lgw-item" title="Stars">
          <span class="lgw-icon">${ICONS.star}</span><span class="lgw-text">${stars}</span>
        </span>
        <span class="lgw-sep"></span>
        <span class="lgw-item" title="Forks">
          <span class="lgw-icon">${ICONS.fork}</span><span class="lgw-text">${forks}</span>
        </span>
      </div>
    `;

        if (isMenuItem) {
            // Wrap in <li> for menu insertion
            const li = document.createElement("li");
            const a = document.createElement("a");
            a.className = "lethe-github-widget";
            a.href = REPO_URL;
            a.target = "_blank";
            a.rel = "noopener noreferrer";
            a.title = `${repoName} on GitHub`;
            a.innerHTML = widgetHTML;
            li.appendChild(a);
            li.id = "lethe-github-widget-container";
            return li;
        } else {
            // Direct anchor element
            const root = document.createElement("a");
            root.id = "lethe-github-widget";
            root.className = "lethe-github-widget";
            root.href = REPO_URL;
            root.target = "_blank";
            root.rel = "noopener noreferrer";
            root.title = `${repoName} on GitHub`;
            root.innerHTML = widgetHTML;
            return root;
        }
    }

    /**
     * Insert widget into the header, or update it if it already exists.
     */
    function insertOrUpdate(data) {
        const { container, before, type } = findInsertTarget();
        if (!container) return;

        const isMenuItem = (type === 'menu');
        const existingId = isMenuItem ? 'lethe-github-widget-container' : 'lethe-github-widget';
        const existing = document.getElementById(existingId);

        // If widget exists and is in the correct container, just update content
        if (existing && existing.parentNode === container) {
            // Check if position is still correct
            const shouldBeBefore = before && before.parentNode === container;
            const isCorrectlyPositioned = shouldBeBefore
                ? existing.nextSibling === before
                : true;

            if (isCorrectlyPositioned) {
                // Update content without DOM manipulation
                const newNode = buildWidget(data, isMenuItem);
                if (isMenuItem) {
                    existing.innerHTML = newNode.innerHTML;
                } else {
                    existing.innerHTML = newNode.innerHTML;
                }
                return;
            }
        }

        // Build new widget
        const node = buildWidget(data, isMenuItem);

        // Remove old widget if it exists anywhere on the page
        const oldWidget = document.getElementById('lethe-github-widget');
        const oldContainer = document.getElementById('lethe-github-widget-container');
        if (oldWidget) oldWidget.remove();
        if (oldContainer) oldContainer.remove();

        // Insert at correct position
        if (before && before.parentNode === container) {
            container.insertBefore(node, before);
        } else {
            // For menu, prepend to show at the start
            if (isMenuItem && container.firstChild) {
                container.insertBefore(node, container.firstChild);
            } else {
                container.appendChild(node);
            }
        }
    }

    // ---------------------------------------------------------------------------
    // Network
    // ---------------------------------------------------------------------------

    /**
     * Fetch JSON from GitHub API.
     */
    async function fetchJson(url) {
        const res = await fetch(url, {
            headers: { "Accept": "application/vnd.github+json" }
        });
        if (!res.ok) throw new Error(`${url} -> ${res.status}`);
        return res.json();
    }

    /**
     * Refresh widget:
     *  1) render cached values immediately (or placeholders),
     *  2) then fetch live values and update + cache them.
     */
    async function refresh() {
        // Paint quickly: cached values or placeholders
        insertOrUpdate(loadCache() || {});

        try {
            const repo = await fetchJson(API_REPO);

            // Keep only what we need
            const data = {
                stargazers_count: repo?.stargazers_count,
                forks_count: repo?.forks_count
            };

            saveCache(data);
            insertOrUpdate(data);
        } catch {
            // If offline / blocked / rate-limited: keep cached or placeholder values
        }
    }

    // ---------------------------------------------------------------------------
    // Startup
    // ---------------------------------------------------------------------------
    function start() {
        refresh();

        // Debounce window resize events to avoid excessive DOM manipulation
        let resizeTimeout;
        window.addEventListener("resize", () => {
            clearTimeout(resizeTimeout);
            resizeTimeout = setTimeout(() => {
                // Re-insert from cache, ensuring correct positioning
                const cachedData = loadCache() || {};
                insertOrUpdate(cachedData);
            }, 100); // Wait 100ms after resize stops before updating
        });
    }

    // Run when DOM is ready
    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", start);
    } else {
        start();
    }
})();