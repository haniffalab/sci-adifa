(function($) {
    $.fn.matrixplot = function(options) {
        var defaults = {
            backgroundColor: "white",                // the canvas background color
            containerId: "canvas-container",         // the data point transparency
            canvasId: "canvas_plot"                  // the color scale usage
        };
        var settings = $.extend({}, defaults, options);
        if (this.length > 1) {
            this.each(function() { $(this).cellatlas(options) });
            return this;
        }
        // private variables
        var width = $('#' + settings.containerId).parent().width();
        var height = $(window).height() - 106;
        var datasetId = this.attr('data-datasetId');
        var xhrPool = [];
        var active = {
            dataset: {},
            bounds: {},
            samples: {},
        };
        var loaded = {
            genes: {}
        };  
        var resources = {};        

        // private methods
        var startLoader = function(id) {
            $("#canvas-loader").html('<div class="btn-group mb-3"><a class="btn btn-white">Loading...</a></div>')
            $("#loader").html('<div class="spinner"><div class="double-bounce1"></div><div class="double-bounce2"></div></div>')
            $("#canvas-controls").hide()
            $("#continuous-legend").empty();
        }
        var endLoader = function(id) {
            $("#canvas-loader").removeClass().empty()
            $("#loader").removeClass().empty()
            $("#canvas-controls").show()
        }
        var showError = function(id) {
            $("#canvas-loader").html('<div class="btn-group mb-3"><a class="btn btn-white">Error</a></div>')
            $("#canvas-controls").hide()
            $("#loader").removeClass().empty()
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
        var render = function() {   
            // load cookies  
            var obsmKey = Cookies.get('ds' + datasetId + '-obsm-key');
            var colorScaleKey = Cookies.get('ds' + datasetId + '-obs-name');
            var colorScaleType = Cookies.get('ds' + datasetId + '-obs-type');

            //Read the data
            var varList = (typeof Cookies.get('ds' + datasetId + '-var-list') === 'undefined') ? [] : JSON.parse(Cookies.get('ds' + datasetId + '-var-list'));
            if(jQuery.isArray(varList) !== -1) {
                markers = varList
            } else {
                markers = ['ALB','AFP','C3','HP','SAA1','RARRES2','LRP1','NR1H4','NNMT','HPD','CES2','C1R','AOX1','GLUL'];
            }

            var myArrayQry = markers.map(function(el, idx) {
                return 'var_names=' + encodeURIComponent(el);
            }).join('&');

            $.when(
                doAjax(API_SERVER + "api/v1/datasets/" + datasetId + "/plotting/matrixplot?groupby=" + colorScaleKey + "&" + myArrayQry).then(function(data) {

                    // Labels of row and columns
                    var myGroups = []
                    var myVars = data.var_names

                    data.categories.forEach(function (item, index) {
                        if ($('#obs-list-' + colorScaleKey.replace(/[^a-zA-Z0-9]/g, '').toLowerCase() + ' input[name="obs-' + item + '"]').is(':checked')) {
                            myGroups.push(item);
                        }
                    });

                    
                    // create svg and set up a y scale, the height value doesn't matter
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

                    var svg = d3.select("#canvas_plot")
                    .append("svg")
                    .attr("width", width + margin.left + margin.right)
                    .attr("height", height + margin.top + margin.bottom)
                    .append("g")
                    .attr("transform",
                            "translate(" + margin.left + "," + margin.top + ")");

                    // Build X scales and axis:
                    var x = d3.scaleBand()
                        .range([ 0, width ])
                        .domain(myVars)
                        .padding(0.02);
                    svg.append("g")
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
                    svg.append("g")
                        .call(d3.axisLeft(y).tickSizeOuter(0));

                    // Build color scale
                    var myColor = d3.scaleSequential()
                        .domain([data.min_value,data.max_value]).interpolator(d3.interpolateViridis);


                        
                    svg.selectAll()
                        .data(Object.entries(data.values_df), function(d) { return d[1] })
                        .enter()
                        .selectAll()
                        .data(function(d, i) { 
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
                        .attr("width", x.bandwidth() )
                        .attr("height", y.bandwidth() )
                        .style("fill", function(d) { return myColor(d[2])} )

                    createLegend(myColor);
                    endLoader();                  
                }, showError));
        }  
        // create continuous color legend
        function createLegend(colorscale) {
            selector_id = '#continuous-legend';

            var legendheight = 200,
                legendwidth = 80,
                margin = {top: 10, right: 60, bottom: 10, left: 0};
        
            var canvas = d3.select(selector_id)
            .style("height", legendheight + "px")
            .style("width", legendwidth + "px")
            .style("position", "absolute")
            .style("top", "10px")
            .style("right", "20px")            
            .append("canvas")
            .attr("height", legendheight - margin.top - margin.bottom)
            .attr("width", 1)
            .style("height", (legendheight - margin.top - margin.bottom) + "px")
            .style("width", (legendwidth - margin.left - margin.right) + "px")
            .style("border", "1px solid #000")
            .style("position", "absolute")
            .style("top", (margin.top) + "px")
            .style("left", (margin.left) + "px")
            .node();
        
            var ctx = canvas.getContext("2d");
        
            var legendscale = d3.scaleLinear()
            .range([legendheight - margin.top - margin.bottom, 1])
            .domain(colorscale.domain());
        
            // image data hackery based on http://bl.ocks.org/mbostock/048d21cf747371b11884f75ad896e5a5
            var image = ctx.createImageData(1, legendheight);
            d3.range(legendheight).forEach(function(i) {
            var c = d3.rgb(colorscale(legendscale.invert(i)));
            image.data[4*i] = c.r;
            image.data[4*i + 1] = c.g;
            image.data[4*i + 2] = c.b;
            image.data[4*i + 3] = 255;
            });
            ctx.putImageData(image, 0, 0);
        
            // A simpler way to do the above, but possibly slower. keep in mind the legend width is stretched because the width attr of the canvas is 1
            // See http://stackoverflow.com/questions/4899799/whats-the-best-way-to-set-a-single-pixel-in-an-html5-canvas
            /*
            d3.range(legendheight).forEach(function(i) {
            ctx.fillStyle = colorscale(legendscale.invert(i));
            ctx.fillRect(0,i,1,1);
            });
            */
        
            var legendaxis = d3.axisRight()
            .scale(legendscale)
            .tickSize(6)
            .ticks(6);
        
            var svg = d3.select(selector_id)
            .append("svg")
            .attr("height", (legendheight) + "px")
            .attr("width", (legendwidth) + "px")
            .style("position", "absolute")
            .style("left", "0px")
            .style("top", "0px")
        
            svg
            .append("g")
            .attr("class", "axis")
            .attr("transform", "translate(" + (legendwidth - margin.left - margin.right - 1) + "," + (margin.top - 1) + ")")
            .call(legendaxis);
        };                
        // // create continuous color legend
        // var createLegend = function(colorscale) {
        //     var margin = {
        //         top: 10,
        //         right: 10,
        //         bottom: 30,
        //         left: 10
        //     };
        //     var legendheight = 50,
        //         legendwidth = $('#canvas-controls').width() + margin.right + margin.left,
        //         selector_id = '#continuous-legend';
                
        //     d3.select(selector_id).selectAll("svg").remove()

        //     var canvas = d3.select(selector_id)
        //         .style("height", legendheight + "px")
        //         .style("width", legendwidth + "px")
        //         .style("position", "absolute")
        //         .style("top", "60px")
        //         .style("right", (30 - margin.left) + "px")
        //         .style("position", "absolute")
        //         .append("canvas")
        //         .attr("height", 1)
        //         .attr("width", legendwidth - margin.left - margin.right)
        //         .style("height", (legendheight - margin.top - margin.bottom) + "px")
        //         .style("width", (legendwidth - margin.left - margin.right) + "px")
        //         .style("border", "1px solid #3d5170")
        //         .style("position", "absolute")
        //         .style("top", (margin.top) + "px")
        //         .style("left", (margin.left) + "px")
        //         .node();

        //     var ctx = canvas.getContext("2d");

        //     var legendscale = d3.scaleLinear()
        //         .range([1, legendwidth - margin.left - margin.right])
        //         .domain(colorscale.domain());

        //     // image data hackery based on http://bl.ocks.org/mbostock/048d21cf747371b11884f75ad896e5a5
        //     var image = ctx.createImageData(legendwidth, 1);
        //     d3.range(legendwidth).forEach(function(i) {
        //         var c = d3.rgb(colorscale(legendscale.invert(i)));
        //         image.data[4 * i] = c.r;
        //         image.data[4 * i + 1] = c.g;
        //         image.data[4 * i + 2] = c.b;
        //         image.data[4 * i + 3] = 255;
        //     });
        //     ctx.putImageData(image, 0, 0);

        //     // A simpler way to do the above, but possibly slower. Keep in mind the legend width is stretched because the width attr of the canvas is 1
        //     // See http://stackoverflow.com/questions/4899799/whats-the-best-way-to-set-a-single-pixel-in-an-html5-canvas
        //     /*
        //     d3.range(legendheight).forEach(function(i) {
        //         ctx.fillStyle = colorscale(legendscale.invert(i));
        //         ctx.fillRect(0,i,1,1);
        //     });
        //     */

        //     var legendaxis = d3.axisBottom()
        //         .scale(legendscale)
        //         .tickSize(6)
        //         .ticks(4);

        //     var svg = d3.select(selector_id)
        //         .append("svg")
        //         .attr("height", (legendheight) + "px")
        //         .attr("width", (legendwidth) + "px")
        //         .style("position", "absolute")
        //         .style("left", "-1px")
        //         .style("top", "-1px")

        //     svg
        //         .append("g")
        //         .attr("class", "axis")
        //         .attr("transform", "translate(" + margin.left + ",20)") // Axis positioning
        //         .call(legendaxis);
        // };
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
                                                .addClass("fa fa-expand")))))
                    .append(
                        $('<div/>')
                            .attr("id", "canvas_plot"))                            
            );
            
            
            startLoader()  
            render()

            return this;
        };
        this.interact = function(el) {

                $(".colourise").removeClass('active');
                $(".btn-gene-select").removeClass("active");

                if (el.id === 'genes'){
                    var colorScaleKey = el.selectedItems[0];
                    var colorScaleId = 0;
                    var colorScaleType = 'gene';                
                } else if ($(el).hasClass('btn-gene-select')) {
                    var colorScaleKey = $(el).text();
                    var colorScaleId = 0;
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
                
                startLoader()  
                render()
        };
        this.genes = function(el) {
            // get cookie data
            var varList = (typeof Cookies.get('ds' + datasetId + '-var-list') === 'undefined') ? [] : JSON.parse(Cookies.get('ds' + datasetId + '-var-list'));
            if(jQuery.isArray(varList) === -1) {
                varList = []
            } 
            if ($(el).hasClass('active')) {
                decolorize();
            } else {
                gene = $(el).text()
                if(jQuery.inArray(gene, varList) === -1) {
                    varList.push(gene)
                } 
                $(el).addClass('active');

                Cookies.set('ds' + datasetId + '-var-list', JSON.stringify(varList), {
                    expires: 30
                })
                
                startLoader();
                render()
            }
        };        
        this.redraw = function() {
            startLoader()  
            render()
        }
        this.removeColor = function() {
            decolorize()
        } 
        this.resetZoom = function() {
            delete currentViewState;
            startLoader()  
            render()
            endLoader()   
        } 
        this.zoomIn = function() {
            const {Viewport} = deck;
            var viewState = (typeof currentViewState != "undefined") ? currentViewState : cellDeck.viewManager.viewState;
            currentViewState = new Viewport({
                width: viewState.width,
                height: viewState.height,
                longitude: viewState.longitude,
                latitude: viewState.latitude,
                zoom: viewState.zoom+1
            });
            cellDeck.setProps({initialViewState: currentViewState});  
        } 
        this.zoomOut = function() {
            const {Viewport} = deck;
            var viewState = (typeof currentViewState != "undefined") ? currentViewState : cellDeck.viewManager.viewState;
            currentViewState = new Viewport({
                width: viewState.width,
                height: viewState.height,
                longitude: viewState.longitude,
                latitude: viewState.latitude,
                zoom: viewState.zoom-1
            });
            cellDeck.setProps({initialViewState: currentViewState});  
        } 
        this.embedding = function(el) {
            Cookies.set('ds' + datasetId + '-obsm-key', $(el).data('name'), {
                expires: 30
            })
            startLoader();
            abort();
            loadData();
        }
        $('.select2-gene-search').select2({
            placeholder: "Search for genes",
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
                .attr("id", "'gene-deg-" + data.id)
                .addClass("btn-gene-select btn btn-outline-info btn-sm")
                    .text(data.id)
            );
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