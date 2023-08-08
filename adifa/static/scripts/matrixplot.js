/* global Cookies */
/* global API_SERVER */
/* global d3 */
(function ($) {
  $.fn.matrixplot = function (options) {
    const defaults = {
      backgroundColor: 'white', // the canvas background color
      containerId: 'matrixplot-container', // the data point transparency
      canvasId: 'canvas_plot' // the color scale usage
    }
    const settings = $.extend({}, defaults, options)
    if (this.length > 1) {
      this.each(function () { $(this).cellatlas(options) })
      return this
    }
    // private variables
    const widthParent = $('#' + settings.containerId).parent().width()
    const datasetId = this.attr('data-datasetId')
    const modality = this.attr('data-modality')
    let xhrPool = []
    const active = {
      dataset: {},
      plot: {},
      dataframe: {}
    }
    let colorScaleKey
    let colorScaleId
    let colorScale
    let varList
    const resources = {}
    // attach a function to be executed before an Ajax request is sent.
    $(document).ajaxSend(function (e, jqXHR, options) {
      xhrPool.push(jqXHR)
    })

    $(document).ajaxComplete(function (e, jqXHR, options) {
      xhrPool = $.grep(xhrPool, function (x) { return x !== jqXHR })
    })

    const abort = function () {
      $.each(xhrPool, function (idx, jqXHR) {
        jqXHR.abort()
      })
    }

    const escapeSelector = function (s) {
      return s.replace(/(:|\.|\[|\])/g, '\\$1')
    }

    // private methods
    const startLoader = function (id) {
      $('#matrixplot-container').show()
      $('#canvas-loader').html('<div class="btn-group mb-3"><a class="btn btn-white">Loading...</a></div>')
      $('#loader').html('<div class="spinner"><div class="double-bounce1"></div><div class="double-bounce2"></div></div>')
      $('#canvas-controls').hide()
      $('#continuous-legend').empty()
      $('#matrixplot-alert').empty().addClass('d-none')
      $('#canvas_plot').empty()
    }

    const endLoader = function (id) {
      $('#canvas-loader').removeClass().empty()
      $('#loader').removeClass().empty()
      $('#canvas-controls').show()
    }

    const showMessage = function (message = null) {
      if (message) {
        $('#matrixplot-alert').removeClass('d-none').append(
          $('<div/>')
            .addClass('alert alert-primary')
            .attr('role', 'alert')
            .text(message))
      }
      $('#canvas-loader').empty()
      $('#canvas-controls').hide()
      $('#matrixplot-container').hide()
      $('#loader').removeClass().empty()
    }

    const showError = function (error) {
      if (!(error.status === 0 && error.statusText === "abort")){
        if (error.statusText) {
          $('#matrixplot-alert').removeClass('d-none').append(
            $('<div/>')
              .addClass('alert alert-danger')
              .attr('role', 'alert')
              .text(error.statusText))
        }
        $('#canvas-loader').html('<div class="btn-group mb-3"><a class="btn btn-white">Error</a></div>')
        $('#canvas-controls').hide()
        $('#loader').removeClass().empty()
      }
    }

    const transform = function (object) {
      return Object
        .entries(object)
        .map(([key, value]) => Object.assign({ key }, value && typeof value === 'object'
          ? { value: '', children: transform(value) }
          : { value, children: [] }
        ))
    }

    const doAjax = function (url, async = true) {
      let result
      // return stored resource
      if (resources[url]) {
        return resources[url]
      }
      // attempt ajax call
      try {
        result = $.ajax({
          url,
          type: 'GET',
          async
        }).done(function () {
          resources[url] = result
        })
          .fail(function (jqXHR, textStatus) {
            return false
          })
        return result
      } catch (error) {
        console.error(error)
        return false
      }
    }

    const loadData = function () {
      if (modality === 'muon') {
        showError('Select a modality with expression matrix to plot.')
        return
      }
      if (!colorScaleKey) {
        showMessage('Please select a group from the list of observations')
        return
      }

      if (!varList.length) {
        showMessage('Select genes of interest from the sidebar on the right')
        return
      }

      const markersQuery = varList.map(function (el, idx) {
        return 'var_names=' + encodeURIComponent(el)
      }).join('&')

      $.when(
        doAjax(API_SERVER + 'api/v1/datasets/' + datasetId + '/plotting/matrixplot?groupby=' + colorScaleKey + '&' + markersQuery + '&modality=' + modality).then(function (data) {
          // update active plot data
          active.data = data
          active.dataframe = transform(data.values_df)
          // get data
          render()
        }, showError))
    }

    const render = function () {
      // ##########################################################################
      // Patrick.Brockmann@lsce.ipsl.fr
      // ##########################################################################

      //= =================================================
      // References
      // http://bl.ocks.org/Soylent/bbff6cc507dca2f48792
      // http://bost.ocks.org/mike/selection/
      // http://bost.ocks.org/mike/join/
      // http://stackoverflow.com/questions/9481497/understanding-how-d3-js-binds-data-to-nodes
      // http://bost.ocks.org/mike/miserables/
      // http://bl.ocks.org/ianyfchang/8119685

      //= =================================================
      $('#matrixplot-container').show()

      const tooltip = d3.select('#canvas_plot')
        .append('div')
        .style('position', 'absolute')
        .style('visibility', 'hidden')

      if (colorScaleKey) {
        $('#color-scale-value').html(colorScaleKey)
        $('#color-scale-value').removeClass('d-none')
        $('#color-scale').removeClass('disabled')
        $('#color-scale-remove').removeClass('d-none')
        if (colorScaleId) {
          $('#collapse' + colorScaleId).collapse('show')
        }
      } else { // decolor
        $('.colourise').removeClass('active')
        $('#color-scale-remove').addClass('d-none')
        $('#color-scale-value').empty().addClass('d-none')
        $('#color-scale').addClass('disabled')
        $('.btn-gene-select').removeClass('active')
      }
      if (varList.length) {
        console.log(active.dataset)
        const varType = active.dataset.modality === 'prot' ? 'protein' : 'genes'
        console.log(varType)
        varList.forEach(function (v) {
          if ($('#gene-deg-' + escapeSelector(v)).length) {
            $('#gene-deg-' + escapeSelector(v)).addClass('active')
          } else {
            $('#search-' + varType + '-selected').append(
              $('<button/>')
                .attr('type', 'button')
                .attr('id', 'gene-deg-' + v)
                .attr('data-gene', v)
                .addClass('btn-gene-select btn btn-outline-info btn-sm active')
                .text(v)
            )
          }
        })
      }

      // Labels of row and columns
      const myVars = active.data.var_names
      let myGroups
      if ($('div[data-name="' + escapeSelector(colorScaleKey) + '"]').data('type') === 'continuous') {
        myGroups = active.data.categories
      } else {
        myGroups = []
        active.data.categories.forEach(function (item, index) {
          if ($('#obs-list-' + colorScaleId + ' input[name="obs-' + escapeSelector(item) + '"]').is(':checked')) {
            myGroups.push(item)
          }
        })
      }

      let svg = d3.select('#canvas_plot')
        .append('svg')

      // Build X scales and axis:
      let y = d3.scaleBand().domain(myGroups)
      let x = d3.scaleBand().domain(myVars)

      // Build X scales and axis:
      const yAxis = svg.append('g')
        .call(d3.axisLeft(y).tickSizeOuter(0))

      // Build X scales and axis:
      const xAxis = svg.append('g')
        .call(d3.axisBottom(x).tickSizeOuter(0))

      xAxis.selectAll('text')
        .style('text-anchor', 'end')
        .attr('dx', '-.8em')
        .attr('dy', '.15em')
        .attr('transform', 'rotate(-65)')

      // determine max width of text label
      let mW = 30
      yAxis.selectAll('.tick').each(function (d) {
        const w = this.getBBox().width
        if (w > mW) mW = w
      })

      let mH = 30
      xAxis.selectAll('.tick').each(function (d) {
        const h = this.getBBox().height
        if (h > mH) mH = h
        // console.log(mH)
      })

      // set the dimensions and margins of the graph
      const margin = { top: 0, right: 0, bottom: mH, left: mW }
      const width = (myVars.length * 20)
      const height = (myGroups.length * 20)

      // append the svg object to the body of the page
      d3.select('#canvas_plot').selectAll('svg').remove()

      // console.log(widthParent)
      svg = d3.select('#canvas_plot')
        .append('svg')
        .style('position', 'relative')
        // .style("width", width + margin.left + margin.right)
        .style('width', '100%')
        .style('display', 'block')
        .style('min-height', height + margin.top + margin.bottom + 130 + 'px')
        .append('g')
        .attr('transform',
          'translate(' + (widthParent / 2 - width / 2) + ',100)')

      // Build X scales and axis:
      x = d3.scaleBand()
        .range([0, width])
        .domain(myVars)
        .padding(0.02)

      svg.append('g')
        .attr('transform', 'translate(0,' + height + ')')
        .call(d3.axisBottom(x).tickSizeOuter(0))
        .selectAll('text')
        .style('text-anchor', 'end')
        .attr('dx', '-.8em')
        .attr('dy', '.15em')
        .attr('transform', 'rotate(-65)')

      // Build X scales and axis:
      y = d3.scaleBand()
        .range([height, 0])
        .domain(myGroups)
        .padding(0.02)

      svg.append('g')
        .call(d3.axisLeft(y).tickSizeOuter(0))

      // $("#canvas_plot").height(height)

      let interpolator
      if (colorScale === 'Turbo') {
        interpolator = d3.interpolateTurbo
      } else if (colorScale === 'Inferno') {
        interpolator = d3.interpolateInferno
      } else if (colorScale === 'Magma') {
        interpolator = d3.interpolateMagma
      } else if (colorScale === 'Plasma') {
        interpolator = d3.interpolatePlasma
      } else if (colorScale === 'Cividis') {
        interpolator = d3.interpolateCividis
      } else if (colorScale === 'Warm') {
        interpolator = d3.interpolateWarm
      } else if (colorScale === 'Cool') {
        interpolator = d3.interpolateCool
      } else {
        interpolator = d3.interpolateViridis
      }

      // Build color scale
      const myColor = d3.scaleSequential()
        .domain([active.data.min_value, active.data.max_value]).interpolator(interpolator)

      let uid = 0
      svg.selectAll()
        .data(Object.entries(active.data.values_df), function (d) { return d[1] })
        .enter()
        .selectAll()
        .data(function (d, i) {
          uid++
          const ownProps = Object.keys(d[1])
          i = ownProps.length
          const resArray = new Array(i) // preallocate the Array
          while (i--) { resArray[i] = [d[0], ownProps[i], d[1][ownProps[i]]] }
          return resArray
        })
        .enter()
        .append('rect')
        .attr('x', function (d) { return x(d[0]) })
        .attr('y', function (d) { return y(d[1]) })
        .attr('class', function (d, i, j) {
          // console.log(uid)
          return 'cell bordered cr' + uid
        })
        .attr('width', x.bandwidth())
        .attr('height', y.bandwidth())
        .on('mouseover', function (event, d) {
          // console.log(event)
          if (d != null) {
            tooltip.html('<div class="heatmap_tooltip">' + d[2] + '</div>')
            tooltip.style('visibility', 'visible')
            tooltip.style('left', (event.offsetX - 44) + 'px')
            tooltip.style('top', (event.offsetY - 50) + 'px')
          } else { tooltip.style('visibility', 'hidden') }
        })
        .on('mouseout', function (event, d) {
          tooltip.style('visibility', 'hidden')
        })

      createLegend(myColor)

      // create svg and set up a y scale, the height value doesn't matter
      const t = svg.transition().duration(500)
      t.selectAll('.cell')
        .style('fill', function (d) { return myColor(d[2]) })

      const zoom = d3.zoom()
        .scaleExtent([0.5, 8])
        .on('zoom', function (event) {
          svg.attr('transform', event.transform)
        })

      svg.call(zoom)

      d3.select('#canvas-zoom-plus')
        .on('click', function () {
          svg.transition().call(zoom.scaleBy, 2)
        })

      d3.select('#canvas-zoom-minus')
        .on('click', function () {
          svg.transition().call(zoom.scaleBy, 0.5)
        })

      d3.select('#canvas-zoom-reset')
        .on('click', function () {
          svg.call(zoom.transform, d3.zoomIdentity.translate((widthParent / 2 - width / 2), 100).scale(1))
        })
      endLoader()
    }

    // create continuous color legend
    function createLegend (colorscale) {
      const selectorId = '#continuous-legend'

      const legendheight = 60
      const legendwidth = 200
      const margin = { top: 0, right: 10, bottom: 40, left: 10 }

      const canvas = d3.select(selectorId)
        .style('height', legendheight + 'px')
        .style('width', legendwidth + 'px')
        .style('top', '20px')
        .style('left', '15px')
        .style('position', 'absolute')
        .append('canvas')
        .attr('height', 1)
        .attr('width', legendwidth - margin.left - margin.right)
        .style('height', (legendheight - margin.top - margin.bottom) + 'px')
        .style('width', (legendwidth - margin.left - margin.right) + 'px')
        .style('border', '1px solid #3d5170')
        .style('position', 'absolute')
        .style('top', (margin.top) + 'px')
        .style('left', (margin.left) + 'px')
        .node()

      const ctx = canvas.getContext('2d')

      const legendscale = d3.scaleLinear()
        .range([1, legendwidth - margin.left - margin.right])
        .domain(colorscale.domain())

      // image data hackery based on http://bl.ocks.org/mbostock/048d21cf747371b11884f75ad896e5a5
      const image = ctx.createImageData(legendwidth, 1)
      d3.range(legendwidth).forEach(function (i) {
        const c = d3.rgb(colorscale(legendscale.invert(i)))
        image.data[4 * i] = c.r
        image.data[4 * i + 1] = c.g
        image.data[4 * i + 2] = c.b
        image.data[4 * i + 3] = 255
      })
      ctx.putImageData(image, 0, 0)

      // A simpler way to do the above, but possibly slower. Keep in mind the legend width is stretched because the width attr of the canvas is 1
      // See http://stackoverflow.com/questions/4899799/whats-the-best-way-to-set-a-single-pixel-in-an-html5-canvas
      /*
          d3.range(legendheight).forEach(function(i) {
              ctx.fillStyle = colorscale(legendscale.invert(i));
              ctx.fillRect(0,i,1,1);
          });
          */

      const legendaxis = d3.axisBottom()
        .scale(legendscale)
        .tickSize(6)
        .ticks(4)

      const svg = d3.select(selectorId)
        .append('svg')
        .attr('height', (legendheight) + 'px')
        .attr('width', (legendwidth) + 'px')
        .style('position', 'absolute')
        .style('left', '-1px')
        .style('top', '-1px')

      svg
        .append('g')
        .attr('class', 'axis')
        .attr('transform', 'translate(' + margin.left + ',20)') // Axis positioning
        .call(legendaxis).append('text')
        .attr('fill', '#3d5170')// set the fill here
        .attr('transform', 'translate(90, 30)')
        .text('Mean expression in group')
    };

    // public methods
    this.initialize = function () {
      // create the explorer tool.
      $('#' + this.attr('id')).append(
        $('<div/>')
          .addClass('card-body')
          .append(
            $('<div/>')
              .attr('id', 'continuous-legend'))
          .append(
            $('<div/>')
              .attr('id', 'loader'))
          .append(
            $('<div/>')
              .attr('id', 'canvas-loader'))
          .append(
            $('<div/>')
              .attr('id', 'canvas-controls')
              .append(
                $('<div/>')
                  .addClass('btn-group mb-3 mr-1')
                  .append(
                    $('<select/>')
                      .attr('id', 'palette')
                      .addClass('form-control')
                      .append(
                        $('<option/>')
                          .attr('value', 'Turbo')
                          .text('Turbo'))
                      .append(
                        $('<option/>')
                          .attr('value', 'Viridis')
                          .text('Viridis'))
                      .append(
                        $('<option/>')
                          .attr('value', 'Inferno')
                          .text('Inferno'))
                      .append(
                        $('<option/>')
                          .attr('value', 'Magma')
                          .text('Magma'))
                      .append(
                        $('<option/>')
                          .attr('value', 'Cividis')
                          .text('Cividis'))
                      .append(
                        $('<option/>')
                          .attr('value', 'Warm')
                          .text('Warm'))
                      .append(
                        $('<option/>')
                          .attr('value', 'Cool')
                          .text('Cool'))))
              .append(
                $('<div/>')
                  .addClass('btn-group mb-3 mr-1')
                  .append(
                    $('<a/>')
                      .attr('id', 'canvas-zoom-plus')
                      .addClass('btn btn-white')
                      .append(
                        $('<i/>')
                          .addClass('fa fa-search-plus')))
                  .append(
                    $('<a/>')
                      .attr('id', 'canvas-zoom-minus')
                      .addClass('btn btn-white')
                      .append(
                        $('<i/>')
                          .addClass('fa fa-search-minus')))
                  .append(
                    $('<a/>')
                      .attr('id', 'canvas-zoom-reset')
                      .addClass('btn btn-white')
                      .append(
                        $('<i/>')
                          .addClass('fa fa-expand')))
                  .append(
                    $('<a/>')
                      .attr('id', 'canvas-gene-reset')
                      .addClass('btn btn-white')
                      .append(
                        $('<i/>')
                          .addClass('fa fa-times')))))
          .append(
            $('<div/>')
              .attr('id', 'canvas_plot'))
      )

      // load cookies
      colorScaleKey = Cookies.get('ds' + datasetId + '-obs-name') || null
      colorScaleId = Cookies.get('ds' + datasetId + '-obs-id') || null
      colorScale = Cookies.get('d3-scale-chromatic') || 'Viridis'
      varList = Cookies.get('ds' + datasetId + '-var-list')
        ? JSON.parse(Cookies.get('ds' + datasetId + '-var-list'))
        : ['ALB', 'AFP', 'C3', 'HP', 'SAA1', 'RARRES2', 'LRP1', 'NR1H4', 'NNMT', 'HPD', 'CES2', 'C1R', 'AOX1', 'GLUL']

      $('#palette').val(colorScale)
      
        // get data
      startLoader()
      $.when(
        doAjax(API_SERVER + 'api/v1/datasets/' + datasetId)).then(function (d) {
        // update active dataset
        active.dataset = d

        if (colorScaleId && !(colorScaleId in active.dataset.data_obs)) {
          colorScaleKey = null
          colorScaleId = null
          Cookies.remove('ds' + datasetId + '-obs-name', {
            path: window.location.pathname
          })
          Cookies.remove('ds' + datasetId + '-obs-id', {
            path: window.location.pathname
          })
        }

        varList = varList.filter(function (n) {
          return active.dataset.data_var.indexOf(n) !== -1
        })
        if (varList.length) {
          Cookies.set('ds' + datasetId + '-var-list', JSON.stringify(varList), {
            expires: 30,
            sameSite: 'Strict',
            path: window.location.pathname
          })
        } else {
          Cookies.remove('ds' + datasetId + '-var-list', {
            path: window.location.pathname
          })
        }

        loadData()
      })

      return this
    }

    this.interact = function (el) {
      $('.colourise').removeClass('active')
      colorScaleKey = $(el).data('name')
      colorScaleId = $(el).data('id')
      Cookies.set('ds' + datasetId + '-obs-name', colorScaleKey, {
        expires: 30,
        sameSite: 'Strict',
        path: window.location.pathname
      })
      Cookies.set('ds' + datasetId + '-obs-id', colorScaleId, {
        expires: 30,
        sameSite: 'Strict',
        path: window.location.pathname
      })

      // get data
      startLoader()
      abort()
      setTimeout(function () { loadData() }, 100)
    }

    this.genes = function (el) {
      // get cookie data
      if (jQuery.isArray(varList) === -1) {
        varList = []
      }
      if ($(el).hasClass('active')) {
        const gene = $(el).attr('data-gene')
        $('#gene-deg-'+escapeSelector(gene)+',#disease-gene-deg-'+escapeSelector(gene)).removeClass('active')
        varList = varList.filter(function (e) { return e !== $(el).data('gene') })
      } else {
        const gene = $(el).text()
        if (jQuery.inArray(gene, varList) === -1) {
          varList.push(gene)
        }
      }
      Cookies.set('ds' + datasetId + '-var-list', JSON.stringify(varList), {
        expires: 30,
        sameSite: 'Strict',
        path: window.location.pathname
      })

      // get data
      startLoader()
      abort()
      setTimeout(function () { loadData() }, 100)
    }

    this.resetGenes = function () {
      $('.btn-gene-select').removeClass('active')
      varList = []
      Cookies.remove('ds' + datasetId + '-var-list', {
        path: window.location.pathname
      })
      // get data
      startLoader()
      loadData()
    }

    this.redraw = function () {
      startLoader()
      render()
      endLoader()
    }

    this.changePalette = function (paletteName) {
      colorScale = paletteName
      Cookies.set('d3-scale-chromatic', colorScale, {
        expires: 30,
        sameSite: 'Strict',
        path: window.location.pathname
      })
      render()
    }

    $('.select2-gene-search').select2({
      placeholder: 'Search by name',
      // closeOnSelect: false,
      sorter: data => data.sort((a, b) => a.text.localeCompare(b.text)),
      ajax: {
        url: API_SERVER + 'api/v1/datasets/' + datasetId + '/search/features',
        data: function (params) {
          const query = {
            search: params.term,
            modality: 'rna',
            type: 'public'
          }

          // Query parameters will be ?search=[term]&type=public
          return query
        }
      }
    }).on('select2:select', function (e) {
      const data = e.params.data
      if (!$('#gene-deg-' + escapeSelector(data.id)).length) {
        $('#search-genes-selected').append(
          $('<button/>')
            .attr('type', 'button')
            .attr('id', 'gene-deg-' + data.id)
            .attr('data-gene', data.id)
            .attr('data-modality', 'rna')
            .addClass('btn-gene-select btn btn-outline-info btn-sm')
            .text(data.id)
        )
        $("button[data-gene='" + escapeSelector(data.id) + "']").trigger('click')
      }
    })

    $('.select2-protein-search').select2({
      // closeOnSelect: false,
      sorter: data => data.sort((a, b) => a.text.localeCompare(b.text)),
      ajax: {
        url: API_SERVER + 'api/v1/datasets/' + datasetId + '/search/features',
        data: function (params) {
          const query = {
            search: params.term,
            modality: 'prot',
            type: 'public'
          }

          // Query parameters will be ?search=[term]&type=public
          return query
        }
      }
    }).on('select2:select', function (e) {
      const data = e.params.data
      if (!$('#gene-deg-' + escapeSelector(data.id)).length) {
        $('#search-protein-selected').append(
          $('<button/>')
            .attr('type', 'button')
            .attr('id', 'gene-deg-' + data.id)
            .attr('data-gene', data.id)
            .attr('data-modality', 'prot')
            .addClass('btn-gene-select btn btn-outline-info btn-sm')
            .text(data.id)
        )
        $("button[data-gene='" + escapeSelector(data.id) + "']").trigger('click')
      }
    })

    $('.select2-disease-search').select2({
      placeholder: 'Search for diseases',
      sorter: data => data.sort((a, b) => a.text.localeCompare(b.text)),
      ajax: {
        url: API_SERVER + 'api/v1/datasets/' + datasetId + '/search/diseases',
        data: function (params) {
          const query = {
            search: params.term,
            type: 'public'
          }

          // Query parameters will be ?search=[term]&type=public
          return query
        }
      }
    }).on('select2:select', function (e) {
      const data = e.params.data
      const genes = data.id.split(',')
      $('#search-genes-disease-set').empty()
      $.each(genes, function (i) {
        let active = false
        if ($('#gene-deg-'+genes[i]).length){
          active = $('#gene-deg-'+genes[i]).hasClass('active')
        }
        $('#search-genes-disease-set').append(
          $('<button/>')
            .attr('type', 'button')
            .attr('id', 'disease-gene-deg-' + genes[i])
            .attr('data-gene', genes[i])
            .addClass('btn-gene-select btn btn-outline-info btn-sm btn-disease-gene')
            .text(genes[i])
            .addClass(active ? 'active' : '')
        )
      })
    })

    return this.initialize()
  }
})(jQuery)
