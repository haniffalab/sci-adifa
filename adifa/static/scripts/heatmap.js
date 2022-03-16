(function($) {
    $.fn.heatmap = function(options) {
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
        var decolorize = function() {
            Cookies.remove('ds' + datasetId + '-obs-name')
            Cookies.remove('ds' + datasetId + '-obs-id')
            Cookies.remove('ds' + datasetId + '-obs-type')
            startLoader()  
            render()
            endLoader() 
        } 
        var abort = function() {
            $.each(xhrPool, function(idx, jqXHR) {
              jqXHR.abort();
            });
          };  
        var getMax = function(arr) {
            let len = arr.length;
            let max = -Infinity;
        
            while (len--) {
                max = arr[len] > max ? arr[len] : max;
            }
            return max;
        };        
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
            var obsmKey = Cookies.get('ds' + datasetId + '-obsm-key');
            var colorScaleKey = Cookies.get('ds' + datasetId + '-obs-name');
            var colorScaleType = Cookies.get('ds' + datasetId + '-obs-type');
            $.when(
                doAjax(API_SERVER + "api/v1/coordinates?embedding=" + obsmKey + "&datasetId=" + datasetId),
                doAjax(API_SERVER + "api/v1/bounds?embedding=" + obsmKey + "&datasetId=" + datasetId)).then(function(a1, a2) {
                    active.samples = a1[0]
                    active.bounds = a2[0] 
                    if (colorScaleKey) {               
                        if (colorScaleType === 'gene'){               
                            var url = API_SERVER + "api/v1/labels?gene=" + colorScaleKey + "&datasetId=" + datasetId;                
                        } else {
                            var url = API_SERVER + "api/v1/labels?obs=" + colorScaleKey + "&datasetId=" + datasetId;                
                        }
                        var requestLabels = doAjax(url, false)
                        if (requestLabels.status === 200) {
                            if (colorScaleType === 'gene'){  
                                active.min = 0;
                                active.max = getMax(requestLabels.responseJSON);
                            }                            
                            active.samples = active.samples.map(function(e, i) {
                                return [e[0], e[1], requestLabels.responseJSON[i]]; 
                            }); 
                        }
                    }
                    render();
                    endLoader();
                    // lazy load obs resources
                    $.each( active.dataset.data_obs, function( key, value ) {
                        doAjax(API_SERVER + "api/v1/labels?obs=" + value.name + "&datasetId=" + datasetId);
                    }); 
                    // Lazy load obsm resources
                    $.each( active.dataset.data_obsm, function( key, value ) {
                        doAjax(API_SERVER + "api/v1/coordinates?embedding=" + value + "&datasetId=" + datasetId);
                        doAjax(API_SERVER + "api/v1/bounds?embedding=" + value + "&datasetId=" + datasetId);                  
                    });                     
                }, showError);
        }
        var render = function() {   
            // init deck
            const {ScatterplotLayer,WebMercatorViewport} = deck;
            // calculate viewport bounding values
            const longitude = (active.bounds['x']['max'] + active.bounds['x']['min']) / 2;
            const latitude = (active.bounds['y']['max'] + active.bounds['y']['min']) / 2;
            const bboxMinLon = Math.max(active.bounds['x']['min'], -179)
            const bboxMinLat = Math.max(active.bounds['y']['min'], -89)
            const bboxMaxLon = Math.min(active.bounds['x']['max'], 179)
            const bboxMaxLat = Math.min(active.bounds['y']['max'], 89)
            // create viewport and fit bounds
            var viewport = new WebMercatorViewport({
                width: width,
                height: height,
                longitude: longitude,
                latitude: latitude,
                zoom: 4
            }).fitBounds([
                [bboxMinLon, bboxMinLat],
                [bboxMaxLon, bboxMaxLat]
            ], {
                padding: {
                    top: 10,
                    bottom: 10,
                    left: 10,
                    right: 10
                }
            });
            // update viewport
            // @TODO: This is causing zoom 0 on init
            // if (typeof currentViewState !== 'undefined') {
            //     var viewport = new WebMercatorViewport({
            //         width: currentViewState.width,
            //         height: currentViewState.height,
            //         longitude: currentViewState.longitude,
            //         latitude: currentViewState.latitude,
            //         zoom: currentViewState.zoom
            //     });
            // }

            currentYear++

            // get cookie data
            var obsmKey = (typeof Cookies.get('ds' + datasetId + '-obsm-key') === 'undefined') ? 'X_umap' : Cookies.get('ds' + datasetId + '-obsm-key');
            var colorScaleKey = (typeof Cookies.get('ds' + datasetId + '-obs-name') === 'undefined') ? null : Cookies.get('ds' + datasetId + '-obs-name');
            var colorScaleId = (typeof Cookies.get('ds' + datasetId + '-obs-id') === 'undefined') ? 0 : Cookies.get('ds' + datasetId + '-obs-id');
            var colorScaleType = (typeof Cookies.get('ds' + datasetId + '-obs-type') === 'undefined') ? null : Cookies.get('ds' + datasetId + '-obs-type');
            // set cell count
            $("#cell-count-value").html(active.samples.length.toLocaleString());
            // set embedding key
            $('#canvas-obsm-key').text(obsmKey);
            if (colorScaleKey) {
                $("#color-scale-value").html(colorScaleKey);
                $("#color-scale-value").removeClass('d-none');
                $("#color-scale").removeClass('disabled');
                $("#color-scale-remove").removeClass('d-none');
            }
            else { // decolor
                $(".colourise").removeClass('active');
                $("#color-scale-remove").addClass('d-none');
                $("#color-scale-value").empty().addClass('d-none');
                $("#color-scale").addClass('disabled');
            }
            if (colorScaleId) {
                $('#collapse' + colorScaleId).collapse("show");
                $('#colourise' + colorScaleId).addClass('active');
            }



            console.log(active.dataset);

            if (colorScaleType == 'categorical') {
                arr = active.dataset['data_obs'][colorScaleKey.replace(/[^a-zA-Z0-9]/g, '').toLowerCase()]['values']
                var catColors = d3.scaleOrdinal().domain($.map(arr, (v, k) => v)).range(["#2f4f4f", "#2e8b57", "#7f0000", "#808000", "#483d8b", "#008000", "#000080", "#8b008b", "#b03060", "#ff0000", "#00ced1", "#ff8c00", "#ffff00", "#00ff00", "#8a2be2", "#00ff7f", "#dc143c", "#00bfff", "#f4a460", "#0000ff", "#f08080", "#adff2f", "#d8bfd8", "#ff00ff", "#1e90ff", "#90ee90", "#ff1493", "#7b68ee", "#ee82ee", "#ffdab9"]);
                var checkboxCheck = new Object();
                for (var k in arr) {
                    if (arr.hasOwnProperty(k)) {
                        //console.log(obs.replace(/[^a-zA-Z0-9]/g, '') + '-' + k)
                        document.getElementById(colorScaleKey.replace(/[^a-zA-Z0-9]/g, '').toLowerCase() + '-' + k).style.backgroundColor = catColors(arr[k]);
                        checkboxCheck[arr[k]] = $('#obs-list-' + colorScaleKey.replace(/[^a-zA-Z0-9]/g, '').toLowerCase() + ' input[name="obs-' + arr[k] + '"]').is(':checked')
                    }
                }
      
                var myColor = function(d) { // Public Method
                    //console.log('#obs-' + d);
                    if (checkboxCheck[d]) {
                        return catColors(d)
                    } else {
                        return "#efefef"
                    }
      
                }
      
                function myRadius(d) {
                    return 100
                    if (checkboxCheck[d]) {
                        return 100
                    } else {
                        return 10
                    }
                }
      
            } else if (colorScaleType == 'continuous') {
                console.log(colorScaleType);
                console.log(active.dataset);
                arr = active.dataset['data_obs'][colorScaleKey.replace(/[^a-zA-Z0-9]/g, '').toLowerCase()]['values']


                var myColor = d3.scaleSequential().domain([active.dataset['data_obs'][colorScaleKey.replace(/[^a-zA-Z0-9]/g, '').toLowerCase()]['min'], active.dataset['data_obs'][colorScaleKey.replace(/[^a-zA-Z0-9]/g, '').toLowerCase()]['max']]).interpolator(d3.interpolateViridis);
                $( "#continuous-legend" ).empty();
                createLegend(myColor);
      
                function myRadius(x) {
                    return 100
                }
            } else if (colorScaleType == 'gene') {
                var myColor = d3.scaleSequential().domain([active.min, active.max]).interpolator(d3.interpolateViridis);
                $( "#continuous-legend" ).empty();
                createLegend(myColor);
      
                function myRadius(x) {
                    return 100
                }
            } else {
                var myColor = function() { // Public Method
                    return '#000';
                }
      
                function myRadius(x) {
                    return 100
                }
            }
            // create data layer
            const layer = new ScatterplotLayer({
              data: active.samples,
              radiusScale: 6,
              radiusMinPixels: 1,
              radiusMaxPixels: 50,
              getPosition: function(d) {
                return [d[0], d[1], 0]
              },
              getFillColor: function(d) {
                v = d3.rgb(myColor(d[2]))
                return [v.r, v.g, v.b]
                return [160, 160, 180, 200]
              },
              getRadius: d => myRadius(d),
              pickable: true, // enable picking (i.e. tooltips)
              updateTriggers: {
                getFillColor: currentYear,
                initialViewState: currentYear,
              }
            })
            // update layer
            cellDeck.setProps({
              layers: [layer] ,
              initialViewState: viewport,
            });             
        }
        // create continuous color legend
        var createLegend = function(colorscale) {
            var margin = {
                top: 10,
                right: 10,
                bottom: 30,
                left: 10
            };
            var legendheight = 50,
                legendwidth = $('#canvas-controls').width() + margin.right + margin.left,
                selector_id = '#continuous-legend';

            var canvas = d3.select(selector_id)
                .style("height", legendheight + "px")
                .style("width", legendwidth + "px")
                .style("position", "absolute")
                .style("top", "60px")
                .style("right", (30 - margin.left) + "px")
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
                .call(legendaxis);
        };
        // public methods 
        this.initialize = function() {
            // load cookies  
            var obsmKey = Cookies.get('ds' + datasetId + '-obsm-key');
            var colorScaleKey = Cookies.get('ds' + datasetId + '-obs-name');
            var colorScaleType = Cookies.get('ds' + datasetId + '-obs-type');

            //Read the data
            d3.json(API_SERVER + "api/v1/datasets/" + datasetId + "/plot/dotplot?obs=" + colorScaleKey + "&markers=ALB,AFP,C3,HP,SAA1,RARRES2,LRP1,NR1H4,NNMT,HPD,CES2,C1R,AOX1,GLUL,CCYP4B1")
            .then(function(data){

                // Labels of row and columns
                var myGroups = data.categories
                var myVars = data.var_names

                // set the dimensions and margins of the graph
                var margin = {top: 30, right: 30, bottom: 100, left: 100},
                width = (myVars.length * 20),
                height = (myGroups.length * 20);

                // append the svg object to the body of the page
                var svg = d3.select("#my_dataviz")
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
                const max = Math.max(Object.values(data.values_df));
                const min = Math.min(Object.values(data.values_df));
                console.log(max)
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
                    .attr("y", function(d) { console.log(d[1]);  return y(d[1]) })
                    .attr("width", x.bandwidth() )
                    .attr("height", y.bandwidth() )
                    .style("fill", function(d) { return myColor(d[2])} )

                createLegend(myColor);

            });
                        
            return this;
        };
        this.colorize = function(el) {
            if ($(el).hasClass('active')) {
                decolorize();
            } else {
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
                
                startLoader();
                abort();
                setTimeout(function(){ loadData(); }, 100); // Defer to improve UX
            }
        };
        this.redraw = function() {
            startLoader()  
            render()
            endLoader()
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