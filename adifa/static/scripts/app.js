(function ($) {
  $(document).ready(function () {
    const escapeSelector = function (s) {
      return s.replace(/(:|\.|\[|\])/g, '\\$1')
    }
    // Enable popovers everywhere.
    $('[data-toggle="popover"]').popover()
    // Enable tooltips everywhere.
    $('[data-toggle="tooltip"]').tooltip()
    // Sidebar toggles
    $('.toggle-sidebar').click(function (e) {
      $('.main-sidebar').toggleClass('open')
    })
    // Init scatterplot
    if ($('#canvas-container').length) {
      const scatterplot = $('#canvas-container').scatterplot()
      // Add events
      // Accordion animations
      $('.obs-values').on('show.bs.collapse', function () {
        $(this).prev('.list-group-item').find('.fa').removeClass('fa-plus-square').addClass('fa-caret-square-up')
        setTimeout(function (el) { scatterplot.colorize(el) }, 100, this) // Defer to improve UX
      }).on('shown.bs.collapse', function () {
        // $('.main-sidebar .nav-wrapper').animate({
        //     scrollTop: $(this).prev(".list-group-item").position().top - 61
        // }, 500, "swing");
      }).on('hide.bs.collapse', function () {
        $(this).prev('.list-group-item').find('.fa').removeClass('fa-caret-square-up').addClass('fa-plus-square')
      })

      $('.colourise').on('click', function (event) {
        if (scatterplot) {
          scatterplot.colorize(this)
          const genes = document.querySelector('#genes')
          if (genes) {
            genes.deselect(genes.selectedItems[0])
          }
        }
      })

      $('#genes').on('update', function (event) {
        if (this.selectedItems[0] && this.selectedItems[0] !== $('#color-scale-value').text() && scatterplot) scatterplot.colorize(this)
      })

      $(document.body).on('click', '.btn-gene-select', function (event) {
        // If selecting from #search-genes-disease-set add to search-genes-selected
        // as selecting another disease will clear #search-genes-disease-set
        const target = $(event.target)
        if ($(target).hasClass("btn-disease-gene")) {
          const gene = target.attr("data-gene")
          if (!$('#search-genes-selected #gene-deg-' + escapeSelector(gene)).length) {
            $('#search-genes-selected').append(
              $('<button/>')
                .attr('type', 'button')
                .attr('id', 'gene-deg-' + gene)
                .attr('data-gene', gene)
                .addClass('btn-gene-select btn btn-outline-info btn-sm')
                .text(gene)
            )
          }
        }
        if ($(this).text() !== $('#color-scale-value').text() && scatterplot) scatterplot.colorize(this)
      })

      $('.obs_value_cb').click(function (event) {
        setTimeout(function () { scatterplot.redraw() }, 100) // Defer to improve UX
      })

      $('.checkall').click(function (event) {
        $('#collapse' + $(this).data('id')).find('input[type=checkbox]').prop('checked', true)
        setTimeout(function () { scatterplot.redraw() }, 100) // Defer to improve UX
      })

      $('.uncheckall').click(function (event) {
        $('#collapse' + $(this).data('id')).find('input[type=checkbox]').prop('checked', false)
        setTimeout(function () { scatterplot.redraw() }, 100) // Defer to improve UX
      })

      $('#color-scale-remove').click(function (event) {
        scatterplot.removeColor(this)
      })

      $('body').on('click', '.canvas-obsm-key', function (event) {
        scatterplot.embedding(this)
      })

      $('#canvas-zoom-reset').click(function (event) {
        scatterplot.resetZoom(this)
      })

      $('#canvas-zoom-plus').click(function (event) {
        scatterplot.zoomIn(this)
      })

      $('#canvas-zoom-minus').click(function (event) {
        scatterplot.zoomOut(this)
      })
    }

    // Init matrixplot
    if ($('#matrixplot-container').length) {
      const matrixplot = $('#matrixplot-container').matrixplot()
      // Add events
      $(document.body).on('change', '#palette', function (event) {
        matrixplot.changePalette($(this).val())
      })
      // Accordion animations
      $('.obs-values').on('show.bs.collapse', function () {
        $(this).prev('.list-group-item').find('.fa').removeClass('fa-plus-square').addClass('fa-caret-square-up')
        setTimeout(function (el) { matrixplot.interact(el) }, 100, this) // Defer to improve UX
      }).on('shown.bs.collapse', function () {
        // $('.main-sidebar .nav-wrapper').animate({
        //     scrollTop: $(this).prev(".list-group-item").position().top - 61
        // }, 500, "swing");
      }).on('hide.bs.collapse', function () {
        $(this).prev('.list-group-item').find('.fa').removeClass('fa-caret-square-up').addClass('fa-plus-square')
      })

      $('.obs_value_cb').click(function (event) {
        setTimeout(function () { matrixplot.redraw() }, 100) // Defer to improve UX
      })

      $('.checkall').click(function (event) {
        $('#collapse' + $(this).data('id')).find('input[type=checkbox]').prop('checked', true)
        setTimeout(function () { matrixplot.redraw() }, 100) // Defer to improve UX
      })

      $('.uncheckall').click(function (event) {
        $('#collapse' + $(this).data('id')).find('input[type=checkbox]').prop('checked', false)
        setTimeout(function () { matrixplot.redraw() }, 100) // Defer to improve UX
      })

      $('#canvas-gene-reset').click(function (event) {
        matrixplot.resetGenes(this)
      })

      $(document.body).on('click', '.btn-gene-select', function (event) {
        // If selecting from #search-genes-disease-set add to search-genes-selected
        // as selecting another disease will clear #search-genes-disease-set
        const target = $(event.target)
        if ($(target).hasClass("btn-disease-gene")) {
          const gene = target.attr("data-gene")
          if (!$('#search-genes-selected #gene-deg-' + escapeSelector(gene)).length) {
            $('#search-genes-selected').append(
              $('<button/>')
                .attr('type', 'button')
                .attr('id', 'gene-deg-' + gene)
                .attr('data-gene', gene)
                .addClass('btn-gene-select btn btn-outline-info btn-sm')
                .text(gene)
            )
          }
        }
        if ($(this).text() !== $('#color-scale-value').text() && matrixplot) matrixplot.genes(this)
      })
    }
  })
})(jQuery)
