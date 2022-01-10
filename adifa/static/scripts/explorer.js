(function($) {
    $.fn.explorer = function(options) {
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
                                                .addClass("fa fa-expand"))))
                            .append(
                                $('<div/>')
                                    .addClass("btn-group mb-3 dropdown")
                                    .append(
                                        $('<button/>')
                                            .attr("id", "canvas-obsm-key")
                                            .attr("type", "button")
                                            .attr("data-toggle", "dropdown")
                                            .attr("aria-haspopup", "true")
                                            .attr("aria-expanded", "false")
                                            .addClass("btn btn-white dropdown-toggle")
                                            .text("Embedding"))
                                    .append(
                                        $('<div/>')
                                            .attr("id", "canvas-obsm-dropdown")
                                            .attr("aria-labelledby", "canvas-obsm-key")
                                            .addClass("dropdown-menu")                                                                                             
                                            .css("z-index","10000"))))                                                                                             
                    .append(
                        $('<div/>')
                            .attr("id", "canvas_plot"))                            
            );
            // set container size
            $('#' + this.attr('id')).parent().height(height)
            $('#' + this.attr('id')).height(height)
            $('#' + this.attr('id')).width(width)
            // load dataset
            startLoader()
            $.when(
                doAjax(API_SERVER + "api/v1/datasets/" + datasetId)).then(function(d1) {
                // update active dataset
                active.dataset = d1    
                // load cookies  
                var colorScaleKey = Cookies.get('ds' + datasetId + '-obs-name')
                var colorScaleType = Cookies.get('ds' + datasetId + '-obs-type') // @TODO 
                // populate embedding options
                $.each(d1.obsm, function(key, obsm) {
                    $('#canvas-obsm-dropdown').append(
                        $('<a/>')
                            .attr("id", "canvas-obsm-key-" + obsm)
                            .attr("href", "#")
                            .attr("data-name", obsm)
                            .addClass("dropdown-item canvas-obsm-key")
                            .text(obsm))
                });                
                // process embedding option
                var obsmKey = (typeof Cookies.get('ds' + datasetId + '-obsm-key') === 'undefined') ? 'X_umap' : Cookies.get('ds' + datasetId + '-obsm-key');
                Cookies.set('ds' + datasetId + '-obsm-key', obsmKey, {
                    expires: 30
                })
                // get data
                loadData();
            }, showError);
            // init deck
            const {DeckGL,WebMercatorViewport} = deck;
            const viewport = new WebMercatorViewport({
                width: width,
                height: height,
                longitude: 0,
                latitude: 0,
                zoom: 0
            })
            deck.log.enable()
            deck.log.level = 3
            currentYear = null // @TODO refactor
            // create deck
            cellDeck = new DeckGL({
                container: 'canvas_plot',
                width: width,
                height: height,
                initialViewState: viewport,
                controller: true,
                getTooltip: ({
                    object
                }) => object && object[2],
                onViewStateChange: ({
                    viewState
                }) => {
                    // we can manipulate the viewState here
                    currentViewState = viewState;
                    console.log(viewState);
                }
            });
            // convert filesize
            const size = $('#ds-size').html()
            var i = Math.floor( Math.log(size) / Math.log(1024) );
            $('#ds-size').html(( size / Math.pow(1024, i) ).toFixed(2) * 1 + ' ' + ['B', 'kB', 'MB', 'GB', 'TB'][i]);

            return this;
        };
        this.colorize = function(el) {
            console.log(el)
            if ($(el).hasClass('active')) {
                decolorize();
            } else {
                $(".colourise").removeClass('active');

                if (el.id === 'genes'){
                    var colorScaleKey = el.selectedItems[0];
                    var colorScaleId = 0;
                    var colorScaleType = 'gene';                
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
        return this.initialize();
    }
})(jQuery);