/* global Cookies */
/* global API_SERVER */
/* global Plotly */
(function ($) {
  $.fn.spatial = function (options) {
    // const defaults = {
    //   backgroundColor: 'white', // the canvas background color
    //   containerId: 'spatial-container' // the data point transparency
    // }
    // const settings = $.extend({}, defaults, options)
    if (this.length > 1) {
      this.each(function () { $(this).cellatlas(options) })
      return this
    }

    const datasetId = this.attr('data-datasetId')
    const spatialModes = ['counts', 'percentage_within_sections', 'percentage_across_sections', 'gene_expression', 'distribution', 'date', 'proportion']
    const colormaps = ['viridis', 'inferno', 'jet', 'RdBu']
    let xhrPool = []

    let colorScaleKey
    let colorScaleType
    let spatialMode
    let prevObsMode
    let colormap

    const active = {
      dataset: {}
    }

    const resources = {}
    // attach a function to be executed before an Ajax request is sent.
    $(document).ajaxSend(function (e, jqXHR, options) {
      xhrPool.push(jqXHR)
    })
    $(document).ajaxComplete(function (e, jqXHR, options) {
      xhrPool = $.grep(xhrPool, function (x) { return x !== jqXHR })
    })

    // const abort = function () {
    //   $.each(xhrPool, function (idx, jqXHR) {
    //     jqXHR.abort()
    //   })
    // }

    // const startLoader = function () {
    //   imgElem.hide()
    //   $('#spatial-loader').show()
    // }

    // const endLoader = function () {
    //   $('#spatial-loader').hide()
    //   imgElem.show()
    // }

    const showError = function () {
      $('#spatial-div').hide()
      $('#spatial-error').show()
      $('#spatial-loader').hide()
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

    this.initialize = function () {
      $('#spatial-div').show()
      $('#spatial-error').hide()
      $('#spatial-loader').hide()

      populateModes()
      populateColormaps()

      $.when(
        doAjax(API_SERVER + 'api/v1/datasets/' + datasetId)).then(function (d1) {
      // update active dataset
        active.dataset = d1

        // load Cookies
        colorScaleKey = Cookies.get('ds' + datasetId + '-obs-name') || null
        colorScaleType = Cookies.get('ds' + datasetId + '-obs-type') || null
        spatialMode = Cookies.get('ds' + datasetId + '-spatial-mode') || null
        colormap = Cookies.get('ds' + datasetId + '-spatial-cm') || colormaps[0]
        if (spatialMode === null && colorScaleType && colorScaleKey) {
          setMode(colorScaleType === 'gene' ? 'gene_expression' : 'counts')
        } else if (spatialMode && colorScaleType && colorScaleKey) {
          setMode(colorScaleType === 'gene' ? 'gene_expression' : spatialMode)
        }
        setColormap(colormap)

        loadPlot()
      }, showError)

      return this
    }

    const populateModes = function () {
      spatialModes.forEach(function (mode) {
        $('#spatial-mode-dropdown').append(
              `<a id="spatial-mode-${mode}" href="#" data-name="${mode}" class="dropdown-item spatial-mode">${mode.replaceAll('_', ' ')}</a>`
        )
      })
      displayModes()
    }

    const displayModes = function () {
      spatialModes.forEach(function (mode) {
        if (colorScaleType === 'gene') {
          $(`#spatial-mode-${mode}`).css('display', (
            mode === 'gene_expression' ? 'block' : 'none'
          ))
        } else if (colorScaleType === 'boolean') {
          $(`#spatial-mode-${mode}`).css('display', (
            mode === 'proportion' ? 'block' : 'none'
          ))
        } else if (colorScaleType === 'date') {
          $(`#spatial-mode-${mode}`).css('display', (
            mode === 'date' ? 'block' : 'none'
          ))
        } else if (colorScaleType === 'continuous') {
          $(`#spatial-mode-${mode}`).css('display', (
            mode === 'distribution' ? 'block' : 'none'
          ))
        } else if (colorScaleType === 'categorical') {
          $(`#spatial-mode-${mode}`).css('display', (
            ['counts', 'percentage_within_sections', 'percentage_across_sections'].includes(mode)
              ? 'block'
              : 'none'
          ))
        } else {
          $(`#spatial-mode-${mode}`).css('display', 'none')
        }
      })
    }

    const displayControls = function () {
      if (spatialMode === 'distribution') {
        $('#spatial-log-scale').show()
      } else {
        $('#spatial-log-scale').hide()
      }
    }

    const populateColormaps = function () {
      colormaps.forEach(function (cm) {
        $('#spatial-colormap-dropdown').append(
              `<a id="spatial-colormap-${cm}" href="#" data-name="${cm}" class="dropdown-item spatial-colormap">${cm}</a>`
        )
      })
    }

    this.redraw = function () {
      loadPlot()
    }

    const loadPlot = function () {
      // startLoader()
      $('#spatial-mode').text(spatialMode ? spatialMode.replaceAll('_', ' ') : '')
      $('#spatial-error').hide()
      $('#spatial-div').show()

      const obsList = []
      if (colorScaleKey && ['categorical', 'boolean', 'date'].includes(colorScaleType)) {
        const arr = active.dataset.data_obs[colorScaleKey.replace(/[^a-zA-Z0-9]/g, '').toLowerCase()].values
        for (const k in arr) {
          if (Object.prototype.hasOwnProperty.call(arr, k)) {
            if ($('#obs-list-' + colorScaleKey.replace(/[^a-zA-Z0-9]/g, '').toLowerCase() + ' input[name="obs-' + arr[k] + '"]').is(':checked')) {
              obsList.push(arr[k])
            }
          }
        }
      }

      const params = []
      if (spatialMode) {
        params.push('mode=' + spatialMode)
        if (spatialMode === 'distribution') {
          if ($('#btn-check-log-scale').prop('checked')) {
            params.push('scale_log=true')
          }
        }
        if (colorScaleKey) {
          if (colorScaleType === 'gene') {
            params.push('plot_value[]=' + [colorScaleKey])
          } else {
            params.push('cat=' + colorScaleKey)
            obsList.forEach(function (o) {
              params.push('plot_value[]=' + o)
            })
          }
        }
      }
      params.push('colormap=' + colormap)

      const paramsStr = params.length ? '?' + params.join('&') : ''

      $.when(
        doAjax(API_SERVER + 'api/v1/datasets/' + datasetId + '/plotting/spatial' + paramsStr
        ).then(function (dataString) {
          // imgElem.attr('src', imgsrc`${data}`)
          const data = JSON.parse(dataString)
          Plotly.newPlot('spatial-plot', data.data, data.layout, { responsive: true, displaylogo: false })
          // endLoader()
        }, showError))
    }

    this.colorize = function (el, active) {
      if (active) {
        this.decolorize()
      } else {
        if (el.id === 'genes') {
          colorScaleKey = el.selectedItems[0]
          colorScaleType = 'gene'
          setMode('gene_expression')
        } else if ($(el).hasClass('btn-gene-select')) {
          colorScaleKey = $(el).text()
          colorScaleType = 'gene'
          setMode('gene_expression')
        } else {
          colorScaleKey = $(el).data('name')
          colorScaleType = $(el).data('type')
          if (colorScaleType === 'boolean') {
            setMode('proportion')
          } else if (colorScaleType === 'date') {
            setMode('date')
          } else if (colorScaleType === 'continuous') {
            setMode('distribution')
          } else {
            setMode(prevObsMode || 'counts')
          }
        }
        displayModes()
        loadPlot()
      }
    }

    this.decolorize = function () {
      colorScaleKey = null
      colorScaleType = null
      spatialMode = null

      Cookies.remove('ds' + datasetId + '-obs-id', {
        path: window.location.pathname
      })
      Cookies.remove('ds' + datasetId + '-obs-name', {
        path: window.location.pathname
      })
      Cookies.remove('ds' + datasetId + '-obs-type', {
        path: window.location.pathname
      })
      Cookies.remove('ds' + datasetId + '-spatial-mode', {
        path: window.location.pathname
      })

      displayModes()
      displayControls()
      loadPlot()
    }

    this.changeMode = function (el) {
      setMode($(el).data('name'))
      loadPlot()
    }

    const setMode = function (mode) {
      if (colorScaleKey) {
        spatialMode = mode
        Cookies.set('ds' + datasetId + '-spatial-mode', spatialMode, {
          expires: 30,
          sameSite: 'Strict',
          path: window.location.pathname
        })
        if (!['gene_expression', 'distribution', 'date', 'proportion'].includes(mode)) {
          prevObsMode = mode
        }
      }
      displayControls()
    }

    this.changeColormap = function (el) {
      setColormap($(el).data('name'))
      loadPlot()
    }

    const setColormap = function (cm) {
      colormap = cm
      colormaps.forEach(function (c) {
        $(`#spatial-colormap-${c}`).removeClass('selected')
      })
      $(`#spatial-colormap-${cm}`).addClass('selected')
      Cookies.set('ds' + datasetId + '-spatial-cm', colormap, {
        expires: 30,
        sameSite: 'Strict',
        path: window.location.pathname
      })
    }

    return this.initialize()
  }
})(jQuery)
