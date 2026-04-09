$(function() {
      // console.log("ping");
      $('[data-toggle="popover"]').popover({
          trigger: 'hover',
              'placement': 'top'
      });

      $.ajaxSetup({
            headers: { 'X-CSRFToken': _csrf_token}
      });


  //// CSRF EXAMPLE FROM DJANGO DOCS
  // function csrfSafeMethod(method) {
  //     // these HTTP methods do not require CSRF protection
  //     return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
  // }
  // $.ajaxSetup({
  //     beforeSend: function(xhr, settings) {
  //         if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
  //             xhr.setRequestHeader("X-CSRFToken", csrftoken);
  //         }
  //     }
  // });

});