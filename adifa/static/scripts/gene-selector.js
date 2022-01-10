export function GeneSelector(idbase) {
    var text = d3.select(`#${idbase}`);
    var multi = d3.select(`#${idbase}-multi`);
    var list = d3.select(`#${idbase}-multi-list`);
    var choices = d3.select(`#${idbase}-multi-choices`);
    var msg = d3.select(`#${idbase}-multi-message`);
    var datalist = d3.select(`#${ text.attr("list") }`)

    var multiple = false;
    var token = 0;
    var did = 0;

    $(choices.node()).on("click", "label", function(e) {
        d3.select(this.parentNode).classed("selected", !d3.select(this.parentNode).classed("selected"))
        $(document).trigger("genelist:updated");
    })

    $(list.node()).on("click", "li", function(e) {
        var item = d3.select(this)
        item.text("").classed("selected", true);
        item.append("label").text(d => d);
        item.append("button").classed("close-icon", true).text("x")
            .on("click", function() { 
                var leaving = d3.select(this.parentNode).remove();
                if (leaving.datum().includes(text.node().value)) {
                    list.append("li")
                        .datum(leaving.datum())
                        .text(d => d);
                    list.selectAll("li").sort();
                }
                $(document).trigger("genelist:updated");
            });

        item.remove();
        choices.append(() => this);
        
        $(document).trigger("genelist:updated");

        if (list.selectAll("li").empty()) {
            listholder.style("display", "none");
        }
    })
    $(`#${idbase}-multi-close`).on("click", function() {
        multi.select(".gsm-list-holder").style("display", "none");
    })

    $(text.node()).on("input", function(ti) {
        window.clearTimeout(token);
        if (did==0) {console.log("not configured"); return;}
        
        if (ti.target.value.length < 3) {
            multi.select(".gsm-list-holder").style("display", "none");
            return;
        }

        this.token = window.setTimeout(function() {
            d3.json(API_SERVER + "api/v1/datasets/" + did + "/genes?term=" + ti.target.value).then((response) => {
                if (response.genes.length == 0) {
                    if (multiple) {
                        msg.text("No genes found.");
                        multi.select(".gsm-list-holder").style("display", "block");
                    } else {
                        datalist.append("option").classed("nodata", true).text("No genes found.").attr("value", "");
                    }
                } else {
                    if (multiple) {
                    response.genes = response.genes.filter(gene => !choices.selectAll("li").data().includes(gene)).sort();

                    list.selectAll("li")
                        .data(response.genes, d => d)
                        .join(
                            enter => enter.append("li")
                                    .text(d => d)
                        )
                        multi.select(".gsm-list-holder").style("display", "block")
                    } else {
                        datalist.selectAll("option")
                            .data(response.genes.sort(), d => d)
                            .join(
                                enter => enter.append("option")
                                    .attr("value", d => d).text(d => d)
                            )
                    }
                }
            })
        }, 500);
    }).on("focus", () => { if (multiple && !list.selectAll("li").empty()) multi.select(".gsm-list-holder").style("display", "block")});

    this.configure = function(d, m) {
        did = d, multiple = m;
        multi.style("display", m?"block":"none");
    }

    function add(gene) {
        var newgene = choices.append("li")
            .classed("selected", true)
            .datum(gene);
            
        newgene.append("label").text(d => d);

        newgene.append("button").classed("close-icon", true).text("x")
            .on("click", function() { 
                var leaving = d3.select(this.parentNode).remove();
                if (leaving.datum().includes(text.node().value)) {
                    list.append("li")
                        .datum(leaving.datum())
                        .text(d => d);
                    list.selectAll("li").sort();
                }
                $(document).trigger("genelist:updated");
            });

    }
    
    this.push = function(genes) {
        if (genes instanceof Array) {
            genes.forEach((gene) => add(gene))
        } else {
            add(genes);
        }

        $(document).trigger("genelist:updated");
    }

    console.log("Gene Selector Created")
}
