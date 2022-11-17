/* global Cookies */
/* global API_SERVER */
(function ($) {
  $.fn.spatial = function (options) {
    const defaults = {
      backgroundColor: 'white', // the canvas background color
      containerId: 'spatial-container', // the data point transparency
    }
    const settings = $.extend({}, defaults, options)
    if (this.length > 1) {
      this.each(function () { $(this).cellatlas(options) })
      return this
    }

    const datasetId = this.attr('data-datasetId')
    const imgElem =  $('#spatial-img')
    let image = ''
    let xhrPool = []

    const imgsrc = function(strings, image){
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

    const abort = function () {
      $.each(xhrPool, function (idx, jqXHR) {
        jqXHR.abort()
      })
    }

    const showError = function (id) {
      $('#canvas-loader').html('<div class="btn-group mb-3"><a class="btn btn-white">Error</a></div>')
      $('#canvas-controls').hide()
      $('#loader').removeClass().empty()
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
      loadPlot();
      return this
    }
    
    const loadPlot = function () {
      $.when(
        doAjax(API_SERVER + 'api/v1/datasets/' + datasetId + '/plotting/spatial').then(function (data) {
          imgElem.attr("src", imgsrc`${data}`)
        }, showError))
    }

    return this.initialize()
  }
})(jQuery)
