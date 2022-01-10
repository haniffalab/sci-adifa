var dotplot = {};
window.addEventListener("DOMContentLoaded", (e) => {

    dotplot.loadtoken = 0;
    dotplot.layout = $('input[name=display-layout]:checked').val();
    dotplot.type = $('input[name=display-tool]:checked').val();
    dotplot.textsize = 16;
    dotplot.boxsize = 20;
    dotplot.colorScale = d3.interpolateRgbBasis(["#98f5ff","#D2691E", "#8b008b"]);
    dotplot.canvas = {};
    dotplot.canvas.width = 800;
    dotplot.canvas.height = 600;
    dotplot.canvas.svg = d3.select("svg#dotplot");
    dotplot.canvas.margin = [40,0,200,220]; //top,right,bottom,left
    dotplot.canvas.base = dotplot.canvas.svg.append("g").attr("transform", "translate(30,30)");
    dotplot.canvas.plot = dotplot.canvas.base.append("g").attr("id", "plot");
    dotplot.canvas.rowlabels = dotplot.canvas.base.append("g").attr("id", "rowlabels");
    dotplot.canvas.collabels = dotplot.canvas.base.append("g").attr("id", "collabels");
    dotplot.dotscale = d3.select("#dphm-size-scale").select("svg");
    dotplot.heatscale = d3.select("#dphm-color-scale").select("svg");
    dotplot.datasetID = dotplot.canvas.svg.node().dataset.did;
    dotplot.cache = {};
    dotplot.data = [];
    dotplot.changes = true;
    dotplot.displaymax = 1;
    dotplot.displaylimit = ISPROTEIN? 15: Infinity;
    dotplot.maxcolumns = Infinity;
    dotplot.hiddencolumns = 0;
    dotplot.startcolumn = 0;
    dotplot.scaledata = [0,0.2,0.4,0.6,0.8,1];
    
    dotplot.dotscale.attr("width", dotplot.boxsize * dotplot.scaledata.length * 1.25)
        .attr("height", dotplot.boxsize + dotplot.textsize + 6)

    dotplot.heatscale.attr("width", dotplot.boxsize * dotplot.scaledata.length  * 1.25)
        .attr("height", dotplot.textsize + 16)

    dotplot.dotscale.selectAll("circle")
        .data(dotplot.scaledata)
        .join(
            enter => enter.append("circle")
        )

    dotplot.heatscale.selectAll("rect")
        .data(dotplot.scaledata)
        .join(
            enter => enter.append("rect")
        )

    dotplot.dotscale.selectAll("text")
        .data(dotplot.scaledata)
        .join(
            enter => enter.append("text").attr("font-family", "Avenir")
            .attr("font-size", dotplot.textsize * 0.75)
            .attr("font-weight", "bold")
            .attr("fill", "black")
        )

    dotplot.heatscale.selectAll("text")
        .data(dotplot.scaledata)
        .join(
            enter => enter.append("text").attr("font-family", "Avenir")
            .attr("font-size", dotplot.textsize * 0.75)
            .attr("font-weight", "bold")
            .attr("fill", "black")
        )
    
    $(".collapse").on('show.bs.collapse', function() {
        $(this).prev(".list-group-item").find(".fa").removeClass("fa-plus-square").addClass("fa-caret-square-up");

    }).on('shown.bs.collapse', function() {
        $('html, body').animate({
            scrollTop: 0
        }, 500, "swing");
        $(document).trigger("category:updated")
    }).on('hide.bs.collapse', function() {
        $(this).prev(".list-group-item").find(".fa").removeClass("fa-caret-square-up").addClass("fa-plus-square");
    });

    $(".obs_value_cb").on("change", function() {
        $(document).trigger("category:updated");
    })

    $(".checkall").click(function(event) {
        $('#collapse' + $(this).data('id')).find('input[type=checkbox]').prop('checked', true);
        $(document).trigger("category:updated");
    });
    
    $(".uncheckall").click(function(event) {
        $('#collapse' + $(this).data('id')).find('input[type=checkbox]').prop('checked', false);
    });

    $('input[name=display-tool]').on("change", function() {
        if (this.checked) {
            dotplot.type = this.value;
            dotplot.refresh();
        }    
    })

    $('input[name=display-layout]').on("change", function() {
        if (this.checked) {
            dotplot.layout = this.value;
            dotplot.init();
            dotplot.refresh();
        }    
    })

    $('.colourise').css("display","none");
    
    window.addEventListener("resize", () => {dotplot.init(); dotplot.refresh(); })

    document.querySelector("#canvas-scroller input[type='range']").addEventListener("input", (e) => {
        dotplot.startcolumn = parseInt(e.target.value);
        dotplot.refresh();
    })

    $(document).on("genelist:updated", function() {
        ensurecategories(); ensuregenes();
        if (dotplot.changes) window.setTimeout(dotplot.load, 100);
    })

    $(document).on("category:updated", function() {
        ensurecategories(); ensuregenes();
        if (dotplot.changes) window.setTimeout(dotplot.load, 100);
    })

    function ensurecategories() {
        var target = $($("#obs-accordion [aria-expanded='true']").attr("data-target")).get(0)
        if (target) {
            if (dotplot.catname == target.dataset.name) {
                var incats = $(target).find("input[type='checkbox']:checked").map((i,e) => e.value).get();
                dotplot.categories = dotplot.categories.filter(c => incats.includes(c)).concat(incats.filter(i => !dotplot.categories.includes(i)));
            } else {
                dotplot.catname =  target.dataset.name;
                dotplot.categories = $(target).find("input[type='checkbox']:checked").map((i,e) => e.value).get();
            }
        } else {
            dotplot.categories = [];
        }
    }

    function ensuregenes() {
        var selectorgenes = document.querySelector("gene-selector").selectedItems;
        if (!dotplot.genes) dotplot.genes = [];
        dotplot.genes = dotplot.genes.filter(g => selectorgenes.includes(g)).concat(selectorgenes.filter(g => !dotplot.genes.includes(g)))
    }

    dotplot.load = function() {

        if (dotplot.categories.length == 0 || dotplot.genes.length == 0) {
            console.log("array lengths: genes: " + dotplot.genes.length + ", categories: " + dotplot.categories.length); return;
        }

        dotplot.canvas.plot.selectAll("*").remove();
        dotplot.displaymax = 1;
        
        dotplot.init();
        dotplot.refresh();

        dotplot.genes.forEach(g => {
            let key = g + ":" + dotplot.catname;

            new Promise((resolve, reject) => {
                if (dotplot.cache[key]) {
                    resolve(dotplot.cache[key])
                } else {
                    var url = API_SERVER.concat("api/v1/datasets/", dotplot.datasetID, "/cxg?obs=", dotplot.catname, "&genes=", g)
                    fetch(url).then(response => response.json(),
                    error => {
                        var geneurl = API_SERVER.concat("api/v1/datasets/", dotplot.datasetID, "/genes?term=", g);
                        fetch(geneurl).then(response => response.json()).then(data => {
                            if (data.genes.includes(g)) {
                                try {
                                    return fetch(url).then(response => response.json());
                                } catch (e) {
                                    console.log(e);
                                }
                            } else {
                                console.log(`Genes "${g} not found in dataset; removing from dotplot`);
                            }
                        })
                    }).then(data => { dotplot.cache[key] = data.data; resolve(data.data)})
                }
            }).then(cached => {
                cached = cached.filter(d => dotplot.categories.includes(d.cat));
                cached.forEach(d => { if (d.expr > dotplot.displaymax) dotplot.displaymax = d.expr });

                dotplot.canvas.plot.selectAll("circle")
                    .data(cached, d => [d.cat, d.gene].join(":"))
                    .join(
                        enter => enter.append("circle")
                            .attr("visibility", "hidden")
/*                          .append("title").text(d => 
`${dotplot.catname}: ${d.cat},
gene: ${d.gene}, 
normalised expression: ${d.expr.toFixed(2)}`
                        )*/,
                        append => {},
                        exit => {}
                    )
    
                dotplot.canvas.plot.selectAll("rect")
                    .data(cached, d => [d.cat, d.gene].join(":"))
                    .join(
                        enter => enter.append("rect").attr("visibility", "hidden")
                          .append("title").text(d => 
`${dotplot.catname}: ${d.cat},
gene: ${d.gene}, 
cells with expression: ${d.expr? (d.count == 0? "< 1": d.count) + "%": "none"},
mean expression: ${d.expr? d.expr.toFixed(2): "none"}`
                        ),
                        append => {},
                        exit => {}
                    )
                dotplot.refresh();
            })
        })
    
    }

    dotplot.init = function() {
        var rows = dotplot.layout[0] == "c"? dotplot.categories: dotplot.genes;
        var cols = dotplot.layout[1] == "c"? dotplot.categories: dotplot.genes;

        dotplot.canvas.rowlabels.selectAll("text")
            .data(rows, d => d)
            .join(
                enter => enter.append("text").text(d => d)
            );

        dotplot.canvas.collabels.selectAll("text")
            .data(cols, d => d)
            .join(
                enter => enter.append("text").text(d => d)
            );
    
        dotplot.canvas.rowlabels.selectAll("text:not(.dragging)")
            .text(d => d)
            .attr("font-family", "Avenir")
            .attr("font-size", dotplot.textsize)
            .attr("fill", "black")
            .call(d3.drag()
                .on("start", function(e,d) {
                    var l = d3.select(this);
                    l.classed("dragging", true)

                    var cumdy = 0, sidx = rows.indexOf(d);

                    e
                        .on("drag", function(ed) { 
                            var dy = ed.dy;

                            cumdy += ed.dy;
                            
                            var cidx = sidx + (Math.trunc((cumdy/dotplot.boxsize)+(dy == 0? 0.2: ((dy/Math.abs(dy)) * 0.2))));
                            cidx = cidx < 0? 0: cidx >= rows.length? rows.length - 1: cidx;
                            
                            if (cidx != rows.indexOf(d)) {
                                rows.splice(cidx, 0, rows.splice(rows.indexOf(d), 1)[0]);
                                dotplot.refresh();
                            }

                            l.attr("y", parseInt(l.attr("y"))+ dy)
                            })
                        .on("end", function() { 
                            l.classed("dragging", false)

                            dotplot.refresh();
                        });
                })
            )

            dotplot.canvas.collabels.selectAll("text:not(.dragging)")
            .text(d => d)
            .attr("font-size", dotplot.textsize)
            .attr("font-family", "Avenir")
            .attr("fill", "black")
            .call(d3.drag()
                .on("start", function(e,d) {
                    var l = d3.select(this);
                    l.classed("dragging", true)
                    var cumdx = 0, sidx = cols.indexOf(d);

                    e
                        .on("drag", function(ed) { 
                            var dx = ed.dx;

                            cumdx += ed.dx;
                            
                            var cidx = sidx + (Math.trunc((cumdx/dotplot.boxsize)+(dx==0? 0.2: (dx/Math.abs(dx)) * 0.2)));
                            cidx = cidx < 0? 0: cidx >= cols.length? cols.length - 1: cidx;
                            
                            if (cidx != cols.indexOf(d)) {
                                cols.splice(cidx, 0, cols.splice(cols.indexOf(d), 1)[0]);
                                dotplot.refresh();
                            }

                            l.attr("x", parseInt(l.attr("x"))+ dx)
                            .attr("transform", function() { 
                                var rotx = parseInt(l.attr("x")) + l.node().getComputedTextLength();
                                return "".concat("rotate(-45,", rotx.toFixed(), ",", 0,")") });
                            })
                        .on("end", function() { 
                            l.classed("dragging", false)

                            dotplot.refresh();
                        });
                })
            )


        // can't set up bottom margin before styling collabels.
        // old version was based on plot height: dotplot.canvas.margin[2] = dotplot.canvas.height - dotplot.boxsize * rows.length;
        // new version measures text, divides by ~root2 (as text is rotated 45 degrees) then adds some padding.
        dotplot.canvas.margin[2] = Math.max(...dotplot.canvas.collabels.selectAll("text").nodes().map(n => n.getComputedTextLength())) / 1.41 + 20;


        // can't set up left margin before styling rowlabels.
        dotplot.canvas.margin[3] = Math.max(...dotplot.canvas.rowlabels.selectAll("text").nodes().map(n => n.getComputedTextLength())) + 10;

        // calculate size of plot v size of svg and container; resize svg appropriately
        //// request to use some kind of scrolling mechanism when plot is too wide for container.

        let svgwidth = dotplot.canvas.margin[1] + dotplot.canvas.margin[3] + (dotplot.boxsize * cols.length) + 30;
        let svgheight = dotplot.canvas.margin[0] + dotplot.canvas.margin[2] + (dotplot.boxsize * rows.length) + 30;

        //// use computedStyles to grab width of dotplot-container and test against svg to see if we need to implement the scrolling feature:
        dotplot.maxcolumns = cols.length;
        let cntwidth = parseFloat(window.getComputedStyle(document.querySelector("#dotplot-container")).width)

        if (svgwidth > cntwidth) {
            let useablewidth = cntwidth - dotplot.canvas.margin[3] - 30; // 30 = x translate of 'base' layer
            
            dotplot.maxcolumns = Math.floor(useablewidth/dotplot.boxsize);
            dotplot.hiddencolumns = cols.length - dotplot.maxcolumns;
            document.querySelector("#canvas-scroller input[type='range']").setAttribute("max", dotplot.hiddencolumns)
            document.querySelector("#canvas-scroller input[type='range']").value = 0
            document.querySelector("#canvas-scroller").style.display = "block";
        } else {
            document.querySelector("#canvas-scroller").style.display = "none";
        }

        // apply svg width and height
        dotplot.canvas.svg
            .attr("height", svgheight)
            .attr("width", svgwidth);

        dotplot.canvas.plot.attr("transform", "".concat("translate(", dotplot.canvas.margin[3], ",", dotplot.canvas.margin[0], ")"))

        dotplot.canvas.rowlabels.attr("transform", "translate(0," + dotplot.canvas.margin[0] + ")");
        dotplot.canvas.collabels.attr("transform", "translate(".concat(dotplot.canvas.margin[3],",",svgheight - dotplot.canvas.margin[2] - 30, ")"))
    }

    dotplot.refresh = function() {
        var rows = dotplot.layout[0] == "c"? dotplot.categories: dotplot.genes;
        var cols = dotplot.layout[1] == "c"? dotplot.categories: dotplot.genes;

        if (dotplot.maxcolumns < cols.length) {
            cols = cols.slice(dotplot.startcolumn, dotplot.maxcolumns + dotplot.startcolumn);
        }

        dotplot.canvas.rowlabels.selectAll("text:not(.dragging)")
            .attr("x", function(d) {return (dotplot.canvas.margin[3] - this.getComputedTextLength() - 5).toFixed() })
            .attr("y", d => (((dotplot.textsize/dotplot.boxsize)+rows.indexOf(d)) * dotplot.boxsize).toFixed())

        dotplot.canvas.collabels.selectAll("text:not(.dragging)")
            .attr("x", function(d) { return ((cols.indexOf(d) * dotplot.boxsize) - this.getComputedTextLength()).toFixed() })
            .attr("y", dotplot.boxsize)
            .attr("transform", function(d) { return "".concat("rotate(-45,", cols.indexOf(d) * dotplot.boxsize, ",", 0,")") })
            .attr("visibility", d => cols.includes(d)? "visible": "hidden")

        dotplot.canvas.plot.selectAll("circle")
            .attr("cx", d => (cols.indexOf(dotplot.layout[1] == "c"? d.cat: d.gene) * dotplot.boxsize) + Math.ceil(dotplot.boxsize/2))
            .attr("cy", d => (rows.indexOf(dotplot.layout[0] == "c"? d.cat: d.gene) * dotplot.boxsize) + Math.ceil(dotplot.boxsize/2))
            .attr("r", d => dotplot.boxsize/2 * d.count/100)
            .attr("fill", d => dotplot.colorScale(d.expr/Math.min(dotplot.displaymax, dotplot.displaylimit)))
            .attr("visibility", d => cols.includes(dotplot.layout[1] == "c"? d.cat: d.gene) ? dotplot.type == "hm"? "hidden":"visible": "hidden")

        dotplot.canvas.plot.selectAll("rect")
            .attr("x", d => (cols.indexOf(dotplot.layout[1] == "c"? d.cat: d.gene) * dotplot.boxsize))
            .attr("y", d => (rows.indexOf(dotplot.layout[0] == "c"? d.cat: d.gene) * dotplot.boxsize))
            .attr("width", dotplot.boxsize)
            .attr("height", dotplot.boxsize)
            .attr("fill", d => dotplot.type == "dp"? "rgba(0,0,0,0)" : d.expr? dotplot.colorScale(d.expr/Math.min(dotplot.displaymax, dotplot.displaylimit)): "#eee")
            //.attr("fill-opacity", d => d.count/100)
            .attr("visibility", d => cols.includes(dotplot.layout[1] == "c"? d.cat: d.gene) ? "visible":"hidden")
                    
        dotplot.dotscale.selectAll("circle")
            .attr("cx", d => (dotplot.scaledata.indexOf(d) * dotplot.boxsize * 1.25) + Math.ceil(dotplot.boxsize/2))
            .attr("cy", Math.ceil(dotplot.boxsize/2))
            .attr("r", d => dotplot.boxsize/2 * d)
            .attr("fill", "#bbb")

        dotplot.dotscale.selectAll("text")
            .text(d => (d*100) + "%")
            .attr("x", function(d) { return (dotplot.scaledata.indexOf(d) * dotplot.boxsize * 1.25) + Math.floor(dotplot.boxsize/2) - Math.floor(this.getComputedTextLength()/2) })
            .attr("y", dotplot.boxsize + dotplot.textsize)

        dotplot.heatscale.selectAll("rect")
            .attr("x", d => (dotplot.scaledata.indexOf(d) * dotplot.boxsize * 1.25))
            .attr("y", 0)
            .attr("width", dotplot.boxsize * 1.25)
            .attr("height", 10)
            .attr("fill", d => dotplot.colorScale(d))

        dotplot.heatscale.selectAll("text")
            .text(d => (d*Math.min(dotplot.displaymax,dotplot.displaylimit)).toFixed(1))
            .attr("x", function(d) { return (dotplot.scaledata.indexOf(d) * dotplot.boxsize * 1.25) + Math.floor(dotplot.boxsize/2) - Math.floor(this.getComputedTextLength()/2) })
            .attr("y", 10 + dotplot.textsize)
            .filter(d => d==1 && dotplot.displaymax > dotplot.displaylimit)
                .text(dotplot.displaylimit.toFixed(1) + "-" + dotplot.displaymax.toFixed(1))

        if (dotplot.displaymax > dotplot.displaylimit) {
            dotplot.heatscale.attr("width", (dotplot.boxsize * dotplot.scaledata.length  * 1.25) + 20)
            dotplot.dotscale.attr("width", (dotplot.boxsize * dotplot.scaledata.length * 1.25) + 20)
        } else {
            dotplot.heatscale.attr("width", (dotplot.boxsize * dotplot.scaledata.length  * 1.25) + 20)
            dotplot.dotscale.attr("width", dotplot.boxsize * dotplot.scaledata.length * 1.25)
        }
        document.querySelector("#dphm-scale").style.display = "inline-flex";
    }

    $($("[data-name=cell.labels]").attr("data-target")).collapse('show');
})