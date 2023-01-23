/* global Cookies */
/* global API_SERVER */
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
    const imgElem = $('#spatial-img')
    let xhrPool = []

    let colorScaleId
    let colorScaleKey
    let colorScaleType

    const active = {
      dataset: {}
    }

    const imgsrc = function (strings, image) {
      return `data:image/png;base64,${image}`
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

      $.when(
        doAjax(API_SERVER + 'api/v1/datasets/' + datasetId)).then(function (d1) {
        // update active dataset
        active.dataset = d1

        // load Cookies
        colorScaleId = Cookies.get('ds' + datasetId + '-obs-id') || null
        colorScaleKey = Cookies.get('ds' + datasetId + '-obs-name') || null
        colorScaleType = Cookies.get('ds' + datasetId + '-obs-type') || null

        loadPlot()
      }, showError)

      loadPlot()
      return this
    }

    this.redraw = function () {
      loadPlot()
    }

    const loadPlot = function () {
      // startLoader()
      $('#spatial-error').hide()
      $('#spatial-div').show()

      const obsList = []
      if (colorScaleKey && colorScaleType === 'categorical') {
        const arr = active.dataset.data_obs[colorScaleKey.replace(/[^a-zA-Z0-9]/g, '').toLowerCase()].values
        for (const k in arr) {
          if (Object.prototype.hasOwnProperty.call(arr, k)) {
            if ($('#obs-list-' + colorScaleKey.replace(/[^a-zA-Z0-9]/g, '').toLowerCase() + ' input[name="obs-' + arr[k] + '"]').is(':checked')) {
              obsList.push(arr[k])
            }
          }
        }
      }

      $.when(
        doAjax(API_SERVER + 'api/v1/datasets/' + datasetId + '/plotting/spatial' +
          (colorScaleKey
            ? ('?cat=' + colorScaleKey) +
            (obsList.length
              ? ('&plot_value[]=' + obsList.join('&plot_value[]='))
              : '') +
            (colorScaleType === 'gene'
              ? ('&gene')
              : '')
            : '')
        ).then(function (data) {
          imgElem.attr('src', imgsrc`${data}`)
          // endLoader()
        }, showError))
    }

    this.colorize = function (el, active) {
      if (active) {
        colorScaleKey = null
        colorScaleId = null
        colorScaleType = null
      } else {
        if (el.id === 'genes') {
          colorScaleKey = el.selectedItems[0]
          colorScaleId = 0
          colorScaleType = 'gene'
        } else if ($(el).hasClass('btn-gene-select')) {
          colorScaleKey = $(el).text()
          colorScaleId = 0
          colorScaleType = 'gene'
        } else {
          colorScaleKey = $(el).data('name')
          colorScaleId = $(el).data('id')
          colorScaleType = $(el).data('type')
        }
      }
      console.log(colorScaleId, colorScaleKey, colorScaleType)
      loadPlot()
    }

    this.decolorize = function () {
      colorScaleId = null
      colorScaleKey = null
      colorScaleType = null
      console.log(colorScaleId, colorScaleKey, colorScaleType)
    }

    return this.initialize()
  }
})(jQuery)
