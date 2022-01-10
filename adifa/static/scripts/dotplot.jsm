

export class DotplotHeatmap extends HTMLElement {

    static get observedAttributes() {
        return ["width", "height", "margins"];
    }

    constructor() {
        super();

        const $ = this.$ = {};
        const svgNS = this.$.svgNS = "http://www.w3.org/2000/svg";

        const root = this.$.root = this.attachShadow({mode: 'open'});

        root.insertAdjacentHTML('beforeend', `          
    <div id="canvas-controls">
        <div class="btn-group btn-group-toggle mb-3" data-toggle="buttons">
          <label class="btn btn-white active">
            <input type="radio" id="display-dotplot" name="display-tool" value="dp" checked>dotplot</label>
          <label class="btn btn-white">
            <input type="radio" id="display-heatmap" name="display-tool" value="hm">heatmap</label>
        </div>               
        <div class="btn-group btn-group-toggle mb-3" data-toggle="buttons">
          <label class="btn btn-white active">
            <input type="radio" id="layout-cellxgene" name="display-layout" value="cg" checked>cells x genes</label>
          <label class="btn btn-white">
            <input type="radio" id="layout-genexcell" name="display-layout" value="gc">genes x cells</label>
        </div>               
    </div>`);

        $.layout = "cg";
        root.querySelectorAll("input[name=display-layout]").forEach(i => { 
            i.addEventListener("change", e => {
                if (e.target.checked) {
                    $.layout = e.target.value;
                    this.refresh();
                }
            })
        });

        this.$.type = "dp";
        root.querySelectorAll("input[name=display-tool]").forEach(i => {
            i.addEventListener("change", e => {
                if (e.target.checked) {
                    $.type = e.target.value;
                    this.refresh();
                }
            })
        });

        $.genes = [];
        $.catn = "";
        $.cats = [];
        $.datacache = {};

        $.svg = document.createElementNS(svgNS, "svg");
        root.appendChild($.svg);
        $.plot = $.svg.appendChild(document.createElementNS(svgNS, "g"));
        $.plot.setAttribute("id", "plot");

        $.rowlabels = $.svg.appendChild(document.createElementNS(svgNS, "g"));
        $.rowlabels.setAttribute("id", "rowlabels");

        $.collabels = $.svg.appendChild(document.createElementNS(svgNS, "g"));
        $.collabels.setAttribute("id", "collabels");

        if (!this.hasAttribute("width")) {
            $.svg.setAttribute("width", 800);
        }

        if (!this.hasAttribute("height")) {
            $.svg.setAttribute("height", 600);
        }

        if (!this.hasAttribute("margins")) {
            this.setAttribute("margins", "80 20 20")
        }

        // mouse event management:
        this.pointer = {};
        this.pointer.down = false;
        this.pointer.target = null;

        this.addEventListener("pointerdown", e => { this.pointer.down == true; this.pointer.target = e.target; })
        this.addEventListener("pointerup", () => { this.pointer.down == false; this.pointer.target = null; })

        this.addEventListener("pointermove", e => {
            if (!this.pointer.down) return;


        })
    }

    attributeChangedCallback(name, oldValue, newValue) {
        switch(name) {
            case "margins":
                let margins = newValue.split(" ").map(v => parseInt(v));

                if (margins.some(v => isNaN(v))) {
                    console.error(`Invalid margins specification "${newValue}" - at least one value is NaN.`)
                    this.$.margins = [0,0,0,0];
                    return;
                }

                if (margins.length > 4) {
                    console.error(`Invalid margins specification "${newValue}" - only 4 values expected.`)
                    this.$.margins = [0,0,0,0];
                    return;
                }

                if (margins.length == 1) {
                    this.$.margins = Array(4).fill(margins[0]);
                } else if (margins.length == 2) {
                    this.$.margins = margins.concat(margins);
                } else if (margins.length == 3) {
                    this.$.margins = margins.concat(margins[1]);
                } else {
                    this.$.margins = margins;
                }
                break;
            case "width":
                this.$.svg.setAttribute("width", newValue)
                break;
            case "height":
                this.$.svg.setAttribute("height", newValue)
                break;

        }
    }

    set dataLoader(value) { this.$.dloader = value; }
    set categoryLoader(value) { this.$.cloader = value; }
    set geneLoader(value) { this.$.gloader = value; }

    get genes() { return this.$.genes }
    set genes(value) { 
        if (Array.isArray(value)) {
            this.$.genes = value;
        }
    }

    addGene(gene, refresh = true) { if (!this.$.genes.includes(gene)) this.$.genes.push(gene); if (refresh) this.refresh() }
    removeGene(gene, refresh = true) { if (this.$,genes.includes(gene)) this.$.genes.splice(this.$.genes.indexOf(gene), 1); if (refresh) this.refresh() }

    get categoryName() { return this.$.catn; }

    set categoryName(value) { 
        if (value != this.$.catn) {
            this.$.catn = value;
            this.$.cats = [];
        }
    }

    get categoryLabels() { return this.$.cats; }
    set categoryLabels(value) { if (Array.isArray(value)) this.$.cats = value }

    addCategoryLabel(label, refresh = true) { if (!this.$.cats.includes(label)) this.$.cats.push(label); if (refresh) this.refresh(); }
    removeCategoryLabel(label, refresh = true) {if (this.$.cats.includes(label)) this.$.cats.splice(this.$.cats.indexOf(label), 1); if (refresh) this.refresh();}
        
    set data(value) { 
        if (Array.isArray(value)) {
            this.$.data = [];
            this.$.cats = [];
            this.$.genes = [];

            value.forEach(d => {
                this.addGene(d.gene);
                this.addCategoryLabel(d.cat);
                this.$.data.push(d);
            });

            this.$.dataloaded = true;
            this.refresh();
        }
    }

    refresh() {
        if (!this.$.dataloaded) {
            if (!this.$.dloader) {
                console.warn("Attempt to refresh dotplot with no data or dataLoader.");
                return;
            }

            this.$.data = [];

            if (this.$.gloader) {
                this.$.genes = this.$.gloader();
            }

            if (this.$.cloader) {
                let cats = this.$.cloader(this.categoryName);
                if (cats.hasOwnProperty('name')) {
                    this.$.catn = cats.name;
                    this.$.cats = cats.labels;
                } else if (this.categoryName && Array.isArray(cats)) {
                    this.$.cats = cats;
                } else {
                    console.error("Unable to resolve category name and labels from categoryLoader.");
                    return;
                }
            }

            let ensureData = g => {
                if (this.$.datacache[g] && this.$.datacache[g][this.categoryName]) {
                    his.$.datacache[g][this.categoryName].forEach(d => {
                        if (this.categoryLabels.includes(d.cat)) {
                            this.$.data.push(d);
                        }
                    });
                    this._draw();
                    this._repaint();
                } else {
                    new Promise((resolve, reject) => {
                        try {
                            resolve(this.$.dloader(g, this.categoryName));
                        } catch(e) {
                            reject(e);
                        }
                    }).then(data => {
                        this.$.datacache[g][this.categoryName] = data;
                        ensureData(g);
                    })
                }

            }
            this.$.genes.forEach(ensureData)
        }
    }

}