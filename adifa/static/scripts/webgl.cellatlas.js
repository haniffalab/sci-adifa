(function($) {
    $.fn.cellatlas = function(options) {
        var defaults = {
            backgroundColor: "white",                // the canvas background color
            particleSize: 2,                         // the size of individual data points
            alphaValue: 1.0,                         // the data point transparency
            paintFeatureScale: 0                     // the color scale usage
        };
        var settings = $.extend({}, defaults, options);
        if (this.length > 1) {
            this.each(function() { $(this).cellatlas(options) });
            return this;
        }
        // private variables
        var canvas = '';                                          // The WebGL canvas element
        var gl = '';                                              // The WebGL context object
        var buffer_data_array = new Float32Array();
        var dr_coordinates = [];
        var type_indices = [];
        // private methods
        var getContext = function(id) {
            var names = ["webgl", "experimental-webgl", "webkit-3d", "moz-webgl"];
            for (var i = 0; i < names.length; i++) {
                try {
                    var gl = canvas.getContext(names[i])
                } catch (e) {}
                if (gl) break
            }
            var vshader = shadersFromScriptElement(gl, "vertex-shader", gl.VERTEX_SHADER),
                fshader = shadersFromScriptElement(gl, "fragment-shader", gl.FRAGMENT_SHADER)
            program = gl.createProgram();
            gl.attachShader(program, vshader)
            gl.attachShader(program, fshader)
            gl.linkProgram(program)
            gl.useProgram(program)
            gl.program = program
            return gl
        }
        var initContext = function(gl) {
            n = buffer_data_array.length / 5
            var vertexColourBuffer = gl.createBuffer()
            gl.bindBuffer(gl.ARRAY_BUFFER, vertexColourBuffer)
      
            var FSIZE = buffer_data_array.BYTES_PER_ELEMENT;
      
            var a_Position = gl.getAttribLocation(gl.program, "a_Position")
            gl.vertexAttribPointer(a_Position, 2, gl.FLOAT, false, FSIZE * 5, 0)
            gl.enableVertexAttribArray(a_Position)
      
            var a_Color = gl.getAttribLocation(gl.program, "a_Color")
            gl.vertexAttribPointer(a_Color, 3, gl.FLOAT, false, FSIZE * 5, 2 * FSIZE)
            gl.enableVertexAttribArray(a_Color)
      
            u_basePointSize = gl.getUniformLocation(gl.program, "u_basePointSize")
            gl.uniform1f(u_basePointSize, settings.particleSize)
      
            u_Alpha = gl.getUniformLocation(gl.program, "u_Alpha")
            gl.uniform1f(u_Alpha, settings.alphaValue)
      
            u_PaintFeatureScale = gl.getUniformLocation(gl.program, "u_PaintFeatureScale")
            gl.uniform1i(u_PaintFeatureScale, settings.paintFeatureScale)
      
            settings.backgroundColor == "white" ? gl.clearColor(1, 1, 1, 1) : gl.clearColor(0, 0, 0, 1)
      
            gl.disable(gl.DEPTH_TEST)
            gl.enable(gl.BLEND)
            gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA)
      
            gl.clear(gl.COLOR_BUFFER_BIT);
            return gl
        }
        var shadersFromScriptElement = function(gl, ID, type) {
            shaderScript = document.getElementById(ID)
            var str = ""
            var k = shaderScript.firstChild;
            while (k) {
                if (k.nodeType == 3) {
                    str += k.textContent;
                }
                k = k.nextSibling
            }
            var shader = gl.createShader(type)
            gl.shaderSource(shader, str)
            gl.compileShader(shader)
            return shader
        }
        var updateBuffer = function() {    
            var buffer_data = [];
            //var colour_by = $('input[name="options"]:checked').val();
            var colour_by = 'cell_type';
            // first update indices to be used - for this read the category control panel radio buttons
            current_indices = []
            for (i = 0; i < document.getElementById("typesControlPanel").childElementCount; i++) {
                if (document.getElementById("typesControlPanel").children[i].children[0].checked) {
                    radio_type = document.getElementById("typesControlPanel").children[i].children[0].id
                    current_indices = current_indices.concat(type_indices[radio_type])
                }
            }
            // now just populate the buffer_data
            if (colour_by == "gene_expression") {
                current_indices.forEach(function(index, i) {
                    buffer_data.push(dr_coordinates[2 * index])
                    buffer_data.push(dr_coordinates[2 * index + 1])
                    buffer_data.push(gene_expression[index])
                    buffer_data.push(gene_expression[index])
                    buffer_data.push(gene_expression[index])
                })
            } else {
                current_indices.forEach(function(index) {
                    buffer_data.push(dr_coordinates[2 * index])
                    buffer_data.push(dr_coordinates[2 * index + 1])
                    buffer_data.push(category_type_colors[3 * index])
                    buffer_data.push(category_type_colors[3 * index + 1])
                    buffer_data.push(category_type_colors[3 * index + 2])
                })
            }
            buffer_data_array = new Float32Array(buffer_data)
            n = buffer_data_array.length / 5
        }
        var draw = function() {    
            updateBuffer()
            if (settings.backgroundColor == "white") {
                gl.clearColor(1, 1, 1, 1)
            } else {
                gl.clearColor(0, 0, 0, 1)
            }
            gl.clear(gl.COLOR_BUFFER_BIT);
            if (n > 0) {
                gl.bufferData(gl.ARRAY_BUFFER, buffer_data_array, gl.STATIC_DRAW)
                gl.drawArrays(gl.POINTS, 0, n)
            }
        }
        var fetchCoordinates = function() {    
            var dr_name_sel = $("input[name='canvas-dr-key']:checked").val();
            //console.log(dr_name_sel);
            $.ajax({
                method: "GET",
                async: false,
                url: "https://developmentcellatlas.ncl.ac.uk/fetch_public_dr_coordinates.php",
                data: { dr_name: dr_name_sel },
                success: function(data) {
                    //console.log( "success" + data );
                    dr_coordinates = [];
                    $.each(data.split(","), function(index, value) {
                        dr_coordinates.push(parseFloat(value));
                    });
                    return dr_coordinates;
                }})
                .done(function() {
                    //console.log( "second success" );
                })
                .fail(function() {
                    //console.log( "error" );
                })
                .always(function() {
                    //console.log( "finished" );
                });

                return dr_coordinates;
        }
        var fetchCategories = function() {    
            var cat_sel = $("#categorySelectMenu").val();
            var loc = window.location.pathname;
            var designfol = window.location.pathname.substring(1, loc.lastIndexOf("/"));
      
            $.ajax({
                method: "GET",
                async: false,
                url: "https://developmentcellatlas.ncl.ac.uk/fetch_public_category.php",
                data: { 
                    category: cat_sel,
                    mot_de_pass: "clear_to_go",
                    design: designfol
                },
                success: function(data) {
                    //console.log( "success" + data );
                    response = data.split("&&");
                    $("#typesControlPanel").html(response[0]);
                    toggleRadio.checked = true;
                    reponse_colors = response[1].split(",");
                    reponse_indices = response[2].split(";");
                    category_type_colors = []
                    type_indices = [];
                    reponse_colors.forEach(function(val) {
                        category_type_colors.push(parseFloat(val));
                    })
                    reponse_indices.forEach(function(indices) {
                        try {
                            indices = indices.split("->");
                            var indices_name = indices[0],
                                indices_values = indices[1].split(",");
                            indices_array = [];
                            indices_values.forEach(function(val) {
                                indices_array.push(parseInt(val))
                            })
                            type_indices[indices_name] = indices_array;
                        } catch (e) {}
                    })
                    return type_indices;
                }})
                .done(function() {
                    //console.log( "second success" );
                })
                .fail(function() {
                    //console.log( "error" );
                })
                .always(function() {
                    //console.log( "finished" );
                });

                return type_indices;
        }
        // public methods 
        this.initialize = function() {
            // define the canvas element.
            canvas = document.getElementById(this.attr('id'));
            // create the renderer
            gl = getContext()
            gl = initContext(gl)
            // maximise the canvas size in the containing element
            var value = parseInt($('#' + this.attr('id') + '-container').width() - 30)
            canvas.width = value
            canvas.height = value
            gl.viewport(0, 0, canvas.width, canvas.height)
            // fetch coordinate data
            dr_coordinates = fetchCoordinates()
            type_indices = fetchCategories()
            // draw to the canvas
            draw();

            return this;
        };
        this.redraw = function() {
            // fetch coordinate data
            dr_coordinates = fetchCoordinates()            
            // draw to the canvas
            draw();
        };
        this.updateCategories = function() {
            // fetch coordinate data
            dr_coordinates = fetchCoordinates()
            type_indices = fetchCategories()
            // draw to the canvas
            draw();
        };
        this.setBackgroundColor = function() {
            this.checked ? $("#" + canvas.id + "-container").css({backgroundColor: "#000"}) : $("#" + canvas.id + "-container").css({backgroundColor: "#fff"});
            settings.backgroundColor = this.checked ? 'dark' : 'white';
            draw()
          }
        return this.initialize();
    }
})(jQuery);