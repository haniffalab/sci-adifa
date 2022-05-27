(function($) {
    $.fn.matrixplot = function(options) {
        var defaults = {
            backgroundColor: "white",                // the canvas background color
            containerId: "matrixplot-container",         // the data point transparency
            canvasId: "canvas_plot"                  // the color scale usage
        };
        var settings = $.extend({}, defaults, options);
        if (this.length > 1) {
            this.each(function() { $(this).cellatlas(options) });
            return this;
        }
        // private variables
        var widthParent = $('#' + settings.containerId).parent().width();
        var heightParent = $(window).height() - 116;
        var datasetId = this.attr('data-datasetId');
        var xhrPool = [];
        var active = {
            plot: {},
            dataframe: {},
        };
        var resources = {};        
        // attach a function to be executed before an Ajax request is sent.
        $(document).ajaxSend(function(e, jqXHR, options){
            xhrPool.push(jqXHR); 
        });
        $(document).ajaxComplete(function(e, jqXHR, options) {
            xhrPool = $.grep(xhrPool, function(x){return x!=jqXHR});
        });
        // private methods
        var startLoader = function(id) {
            $("#canvas-loader").html('<div class="btn-group mb-3"><a class="btn btn-white">Loading...</a></div>')
            $("#loader").html('<div class="spinner"><div class="double-bounce1"></div><div class="double-bounce2"></div></div>')
            $("#canvas-controls").hide()
            $("#continuous-legend").empty();
            $('#matrixplot-alert').empty().addClass('d-none')
            $('#canvas_plot').empty()
        }
        var endLoader = function(id) {
            $("#canvas-loader").removeClass().empty()
            $("#loader").removeClass().empty()
            $("#canvas-controls").show()
        }
        var showError = function(message=false) {
            if (message)
                $('#matrixplot-alert').removeClass('d-none').append(
                    $('<div/>')
                        .addClass("alert alert-primary")
                        .attr("role", "alert")
                        .text(message))
            $("#canvas-loader").html('<div class="btn-group mb-3"><a class="btn btn-white">Error</a></div>')
            $("#canvas-controls").hide()
            $("#loader").removeClass().empty()
        }
        var abort = function() {
            $.each(xhrPool, function(idx, jqXHR) {
              jqXHR.abort();
            });
          };      
        var transform = function(object) {
            return Object
                .entries(object)
                .map(([key, value]) => Object.assign({ key }, value && typeof value === 'object'
                    ? { value: '', children: transform(value) }
                    : { value, children: [] }
                ));
        }
        var validate = function() {
            // load cookies  
            var varList = (typeof Cookies.get('ds' + datasetId + '-var-list') === 'undefined') ? [] : JSON.parse(Cookies.get('ds' + datasetId + '-var-list'));
            if(jQuery.isArray(varList) !== -1) {
                cMarkers = varList
            } else {
                cMarkers = ['ALB','AFP','C3','HP','SAA1','RARRES2','LRP1','NR1H4','NNMT','HPD','CES2','C1R','AOX1','GLUL'];
            }
            // validate markers
            var markers = cMarkers.filter(function(marker) {
                if ($("button[data-gene='" + marker + "']").length) {
                    $("button[data-gene='" + marker + "']").addClass("active")
                    return true;
                }
                return false;
            });
            return markers
        }
        var doAjax = function(url, async=true) {
            let result;
            // return stored resource
            if (resources[url]) {
                return resources[url];
            }
            // attempt ajax call
            try {
                result = $.ajax({
                    url: url,
                    type: 'GET',
                    async: async
                }).done(function() {
                    resources[url] = result
                })
                .fail(function( jqXHR, textStatus ) {
                    return false;
                });
                return result;
            } catch (error) {
                console.error(error);
                return false;
            }
        }    
        var loadData = function() {
            // load cookies  
            if (typeof Cookies.get('ds' + datasetId + '-obs-name') === 'undefined') {
                showError('Please select a group from the list of observations')
                return
            } 
            // get validated markers
            var markers = validate()
            var markersQuery = markers.map(function(el, idx) {
                return 'var_names=' + encodeURIComponent(el);
            }).join('&');

            if (!markers.length > 0) {
                showError('Select genes of interest from the sidebar on the right')
                return
            } 

            $.when(
                doAjax(API_SERVER + "api/v1/datasets/" + datasetId + "/plotting/matrixplot?groupby=" + Cookies.get('ds' + datasetId + '-obs-name') + "&" + markersQuery).then(function(data) {
                // update active plot data
                active.data = data
                active.dataframe = transform(data.values_df)
                // get data
                render();
            }, showError));

        }
        var render = function() {
            // load cookies  
            var obsmKey = Cookies.get('ds' + datasetId + '-obsm-key');
            var colorScaleKey = Cookies.get('ds' + datasetId + '-obs-name');
            var colorScaleType = Cookies.get('ds' + datasetId + '-obs-type');
            var colorScale = (typeof Cookies.get('d3-scale-chromatic') === 'undefined') ? "Viridis" : Cookies.get('d3-scale-chromatic');

            //##########################################################################
            // Patrick.Brockmann@lsce.ipsl.fr
            //##########################################################################
            
            //==================================================
            // References
            // http://bl.ocks.org/Soylent/bbff6cc507dca2f48792
            // http://bost.ocks.org/mike/selection/
            // http://bost.ocks.org/mike/join/
            // http://stackoverflow.com/questions/9481497/understanding-how-d3-js-binds-data-to-nodes
            // http://bost.ocks.org/mike/miserables/
            // http://bl.ocks.org/ianyfchang/8119685

            //==================================================
            var tooltip = d3.select("#canvas_plot")
            .append("div")
            .style("position", "absolute")
            .style("visibility", "hidden");

            // Labels of row and columns
            var myVars = active.data.var_names
            if ($('div[data-name="' + colorScaleKey + '"]').data('type') == "continuous") {
                var myGroups = active.data.categories
            }
            else {
                var myGroups = []
                active.data.categories.forEach(function (item, index) {
                    if ($('#obs-list-' + colorScaleKey.replace(/[^a-zA-Z0-9]/g, '').toLowerCase() + ' input[name="obs-' + item + '"]').is(':checked')) {
                        myGroups.push(item);
                    }
                });
            }
                    





                    var svg = d3.select("#canvas_plot")
                    .append("svg");



                    // Build X scales and axis:
                    var y = d3.scaleBand().domain(myGroups);
                    var x = d3.scaleBand().domain(myVars);

                    // Build X scales and axis:
                    var yAxis = svg.append("g")
                    .call(d3.axisLeft(y).tickSizeOuter(0));

                    // Build X scales and axis:
                    var xAxis = svg.append("g")
                        .call(d3.axisBottom(x).tickSizeOuter(0))

                    xAxis.selectAll("text")  
                        .style("text-anchor", "end")
                        .attr("dx", "-.8em")
                        .attr("dy", ".15em")
                        .attr("transform", "rotate(-65)")

                    // determine max width of text label
                    var mW = 30;
                    yAxis.selectAll(".tick").each(function(d) {
                        var w = this.getBBox().width;
                        if (w > mW) mW = w;
                    });

                    var mH = 30;
                    xAxis.selectAll(".tick").each(function(d) {
                        var h = this.getBBox().height;
                        if (h > mH) mH = h;
                        console.log(mH)
                    });

                    // set the dimensions and margins of the graph
                    var margin = {top: 0, right: 0, bottom: mH, left: mW},
                    width = (myVars.length * 20),
                    height = (myGroups.length * 20);

                    // append the svg object to the body of the page
                    d3.select("#canvas_plot").selectAll("svg").remove()

                    console.log(widthParent)
                    var svg = d3.select("#canvas_plot")
                    .append("svg")
                        .style("position", "relative")
                        //.style("width", width + margin.left + margin.right)
                        .style("width", "100%")
                        .style("display", "block")
                        .style("min-height", height + margin.top + margin.bottom + 130)
                        .append("g")
                            .attr("transform",
                                    "translate("+ (widthParent/2 - width/2) +",100)");

                    // Build X scales and axis:
                    var x = d3.scaleBand()
                        .range([ 0, width ])
                        .domain(myVars)
                        .padding(0.02);
                    var colLabels = svg.append("g")
                        .attr("transform", "translate(0," + height + ")")
                        .call(d3.axisBottom(x).tickSizeOuter(0))
                        .selectAll("text")  
                        .style("text-anchor", "end")
                        .attr("dx", "-.8em")
                        .attr("dy", ".15em")
                        .attr("transform", "rotate(-65)")

                    // Build X scales and axis:
                    var y = d3.scaleBand()
                        .range([ height, 0 ])
                        .domain(myGroups)
                        .padding(0.02);
                    var rowLabels = svg.append("g")
                        .call(d3.axisLeft(y).tickSizeOuter(0));

                    if (colorScale == "Turbo") {
                        interpolator = d3.interpolateTurbo
                    } else if (colorScale == "Inferno") {
                        interpolator = d3.interpolateInferno
                    } else if (colorScale == "Magma") {
                        interpolator = d3.interpolateMagma
                    } else if (colorScale == "Plasma") {
                        interpolator = d3.interpolatePlasma
                    } else if (colorScale == "Cividis") {
                        interpolator = d3.interpolateCividis
                    } else if (colorScale == "Warm") {
                        interpolator = d3.interpolateWarm
                    } else if (colorScale == "Cool") {
                        interpolator = d3.interpolateCool
                    } else {
                        interpolator = d3.interpolateViridis
                    }

                    // Build color scale
                    var myColor = d3.scaleSequential()
                        .domain([active.data.min_value,active.data.max_value]).interpolator(interpolator);

                    uid = 0
                    svg.selectAll()
                        .data(Object.entries(active.data.values_df), function(d) { return d[1] })
                        .enter()
                        .selectAll()
                        .data(function(d, i) { 
                            uid++;
                            var ownProps = Object.keys( d[1] ),
                            i = ownProps.length,
                            resArray = new Array(i); // preallocate the Array
                            while (i--)
                                resArray[i] = [d[0], ownProps[i], d[1][ownProps[i]]];
                            return resArray })
                        .enter()
                        .append("rect")
                        .attr("x", function(d) { return x(d[0]) })
                        .attr("y", function(d) { return y(d[1]) })
                        .attr("class", function(d, i, j) {
                            console.log(uid)
                            return "cell bordered cr" + uid;
                        })                        
                        .attr("width", x.bandwidth() )
                        .attr("height", y.bandwidth() )
                        .on('mouseover', function(event,d) {
                            console.log(event)
                            if (d != null) {
                                tooltip.html('<div class="heatmap_tooltip">' + d[2] + '</div>');
                                tooltip.style("visibility", "visible");
                                tooltip.style("left", (event.offsetX - 44) + "px");
                                tooltip.style("top", (event.offsetY - 50) + "px");
                            } else
                                tooltip.style("visibility", "hidden");
                        })
                        .on('mouseout', function(event,d) {
                            tooltip.style("visibility", "hidden");
                        });
                        
                        createLegend(myColor);

                    // create svg and set up a y scale, the height value doesn't matter
                    var t = svg.transition().duration(500);
                    t.selectAll(".cell")
                        .style("fill", function(d) { return myColor(d[2])} )


                    var zoom = d3.zoom()
                    .scaleExtent([0.5, 8])
                    .on('zoom', function(event) {
                        svg.attr('transform', event.transform);
                    });
                
                    
                    svg.call(zoom);

                    d3.select('#canvas-zoom-plus')
                    .on('click', function() { 
                        svg.transition().call(zoom.scaleBy, 2)
                     });

                     d3.select('#canvas-zoom-minus')
                     .on('click', function() { 
                         svg.transition().call(zoom.scaleBy, 0.5)
                      });

                      d3.select('#canvas-zoom-reset')
                      .on('click', function() { 
                          svg.transition().call(zoom.scaleBy, 1)
                       });
                    endLoader();                  
        }  
        // create continuous color legend
        function createLegend(colorscale) {
            selector_id = '#continuous-legend';

            var legendheight = 60,
                legendwidth = 200,
                margin = {top: 0, right: 10, bottom: 40, left: 10};
        
            var canvas = d3.select(selector_id)
            .style("height", legendheight + "px")
            .style("width", legendwidth + "px")
            .style("top", "20px")
            .style("left", "15px")
            .style("position", "absolute")
            .append("canvas")
            .attr("height", 1)
            .attr("width", legendwidth - margin.left - margin.right)
            .style("height", (legendheight - margin.top - margin.bottom) + "px")
            .style("width", (legendwidth - margin.left - margin.right) + "px")
            .style("border", "1px solid #3d5170")
            .style("position", "absolute")
            .style("top", (margin.top) + "px")
            .style("left", (margin.left) + "px")
            .node();

            var ctx = canvas.getContext("2d");

            var legendscale = d3.scaleLinear()
                .range([1, legendwidth - margin.left - margin.right])
                .domain(colorscale.domain());

            // image data hackery based on http://bl.ocks.org/mbostock/048d21cf747371b11884f75ad896e5a5
            var image = ctx.createImageData(legendwidth, 1);
            d3.range(legendwidth).forEach(function(i) {
                var c = d3.rgb(colorscale(legendscale.invert(i)));
                image.data[4 * i] = c.r;
                image.data[4 * i + 1] = c.g;
                image.data[4 * i + 2] = c.b;
                image.data[4 * i + 3] = 255;
            });
            ctx.putImageData(image, 0, 0);

            // A simpler way to do the above, but possibly slower. Keep in mind the legend width is stretched because the width attr of the canvas is 1
            // See http://stackoverflow.com/questions/4899799/whats-the-best-way-to-set-a-single-pixel-in-an-html5-canvas
            /*
            d3.range(legendheight).forEach(function(i) {
                ctx.fillStyle = colorscale(legendscale.invert(i));
                ctx.fillRect(0,i,1,1);
            });
            */

            var legendaxis = d3.axisBottom()
                .scale(legendscale)
                .tickSize(6)
                .ticks(4);

            var svg = d3.select(selector_id)
                .append("svg")
                .attr("height", (legendheight) + "px")
                .attr("width", (legendwidth) + "px")
                .style("position", "absolute")
                .style("left", "-1px")
                .style("top", "-1px")

            svg
                .append("g")
                .attr("class", "axis")
                .attr("transform", "translate(" + margin.left + ",20)") // Axis positioning
                .call(legendaxis).append("text")
                .attr("fill", "#3d5170")//set the fill here
                .attr("transform","translate(90, 30)")
                .text("Mean expression in group");;
        };                
        // public methods 
        this.initialize = function() {
            // create the explorer tool.
            $('#' + this.attr('id')).append(
                $('<div/>')
                    .addClass("card-body")
                    .append(
                        $('<div/>')
                            .attr("id", "continuous-legend"))
                    .append(
                        $('<div/>')
                            .attr("id", "loader"))
                    .append(
                        $('<div/>')
                            .attr("id", "canvas-loader"))
                    .append(
                        $('<div/>')
                            .attr("id", "canvas-controls")    
                            .append(
                                $('<div/>')
                                    .addClass("btn-group mb-3 mr-1")
                                    .append(
                                        $('<select/>')
                                            .attr("id", "palette")
                                            .addClass("form-control")
                                            .append(
                                                $('<option/>')
                                                .attr("value", "Turbo")
                                                .text("Turbo"))
                                            .append(
                                                $('<option/>')
                                                .attr("value", "Viridis")
                                                .text("Viridis"))
                                            .append(
                                                $('<option/>')
                                                .attr("value", "Inferno")
                                                .text("Inferno"))
                                            .append(
                                                $('<option/>')
                                                .attr("value", "Magma")
                                                .text("Magma"))
                                            .append(
                                                $('<option/>')
                                                .attr("value", "Cividis")
                                                .text("Cividis"))
                                            .append(
                                                $('<option/>')
                                                .attr("value", "Warm")
                                                .text("Warm"))
                                            .append(
                                                $('<option/>')
                                                .attr("value", "Cool")
                                                .text("Cool"))))     
                            .append(
                                $('<div/>')
                                    .addClass("btn-group mb-3 mr-1")
                                    .append(
                                        $('<a/>')
                                            .attr("id", "canvas-zoom-plus")
                                            .addClass("btn btn-white")
                                            .append(
                                                $('<i/>')
                                                .addClass("fa fa-search-plus")))
                                    .append(
                                        $('<a/>')
                                            .attr("id", "canvas-zoom-minus")
                                            .addClass("btn btn-white")
                                            .append(
                                                $('<i/>')
                                                .addClass("fa fa-search-minus")))
                                    .append(
                                        $('<a/>')
                                            .attr("id", "canvas-zoom-reset")
                                            .addClass("btn btn-white")
                                            .append(
                                                $('<i/>')
                                                .addClass("fa fa-expand")))                                                
                                    .append(
                                        $('<a/>')
                                            .attr("id", "canvas-gene-reset")
                                            .addClass("btn btn-white")
                                            .append(
                                                $('<i/>')
                                                .addClass("fa fa-times")))))
                    .append(
                        $('<div/>')
                            .attr("id", "canvas_plot"))
            );
            // get data
            startLoader()
            loadData();
            return this;
        };
        this.interact = function(el) {
            $(".colourise").removeClass('active');
            $(".btn-gene-select").removeClass("active");

            if (el.id === 'genes'){
                var colorScaleKey = el.selectedItems[0];
                var colorScaleId = 0;
                var colorScaleType = 'gene';                
                    var colorScaleType = 'gene';                
                var colorScaleType = 'gene';                
            } else if ($(el).hasClass('btn-gene-select')) {
                var colorScaleKey = $(el).text();
                var colorScaleId = 0;
                var colorScaleType = 'gene';  
                    var colorScaleType = 'gene';  
                var colorScaleType = 'gene';  
                $(el).addClass('active');
            } else {
                var colorScaleKey = $(el).data('name');
                var colorScaleId = $(el).data('id');
                var colorScaleType = $(el).data('type');
                $('#colourise' + colorScaleId).addClass('active');
            }
            Cookies.set('ds' + datasetId + '-obs-name', colorScaleKey, {
                expires: 30
            })
            Cookies.set('ds' + datasetId + '-obs-id', colorScaleId, {
                expires: 30            
            })        
            Cookies.set('ds' + datasetId + '-obs-type', colorScaleType, {
                expires: 30
            })
             
            // get data
            startLoader()
                startLoader()  
            startLoader()
            loadData();
        };
        this.genes = function(el) {
            // get cookie data
            var varList = (typeof Cookies.get('ds' + datasetId + '-var-list') === 'undefined') ? [] : JSON.parse(Cookies.get('ds' + datasetId + '-var-list'));
            if(jQuery.isArray(varList) === -1) {
                varList = []
            } 
            if ($(el).hasClass('active')) {
                $(el).removeClass('active');
                var varFiltered = varList.filter(function(e, y) { return e !== $(el).data('gene') })
                varList = varFiltered
            } else {
                gene = $(el).text()
                if(jQuery.inArray(gene, varList) === -1) {
                    varList.push(gene)
                } 
                $(el).addClass('active');
            }
            Cookies.set('ds' + datasetId + '-var-list', JSON.stringify(varList), {
                expires: 30
            })
            
            // get data
            startLoader()
            loadData();
        };        
        this.resetGenes = function() {
            $(".btn-gene-select").removeClass('active');
            Cookies.remove('ds' + datasetId + '-var-list')
            // get data
            startLoader()
            loadData(); 
        } 
        this.redraw = function() {
            startLoader()  
            render()
            endLoader()
        }
        this.changePalette = function(paletteName) {
            Cookies.set('d3-scale-chromatic', paletteName, {
                expires: 30
            })
            render()
        }
        $('.select2-gene-search').select2({
            placeholder: "Search by name",
            closeOnSelect: false,
            sorter: data => data.sort((a, b) => a.text.localeCompare(b.text)),
            ajax: {
                url: API_SERVER + "api/v1/datasets/" + datasetId + "/search/genes",
                data: function (params) {
                var query = {
                    search: params.term,
                    type: 'public'
                }

                // Query parameters will be ?search=[term]&type=public
                return query;
                }
            }
        }).on('select2:select', function (e) {
            var data = e.params.data;
            $('#search-genes-selected').append(
                $('<button/>')
                .attr("type", "button")
                .attr("id", "gene-deg-" + data.id)
                .attr("data-gene", data.id)
                .addClass("btn-gene-select btn btn-outline-info btn-sm")
                    .text(data.id)
            );
            $("button[data-gene='" + data.id + "']").trigger( "click" );
        });
        $('.select2-disease-search').select2({
            placeholder: "Search for diseases",
            sorter: data => data.sort((a, b) => a.text.localeCompare(b.text)),
            ajax: {
                url:  API_SERVER + "api/v1/datasets/" + datasetId + "/search/diseases",
                data: function (params) {
                var query = {
                    search: params.term,
                    type: 'public'
                }

                // Query parameters will be ?search=[term]&type=public
                return query;
                }
            }
        }).on('select2:select', function (e) {
            var data = e.params.data;
            var genes = data.id.split(",");
            $.each(genes,function(i){
                $('#search-genes-disease-set').append(
                $('<button/>')
                    .attr("type", "button")
                    .attr("id", "'gene-deg-" + genes[i])
                    .addClass("btn-gene-select btn btn-outline-info btn-sm")
                    .text(genes[i])
                );
            });
        });
        return this.initialize();
    }
})(jQuery);