
function showSimilarSamples(div, similar_samples){
    $(div).empty();
    $.each(similar_samples, function(){
        $(div).append('<img src="http://localhost:5000/img_api/svg_fingerprints/' + this + '.svg" /><br/>'); 
    });
}

var query_profile = "";

$(function() {
    // global AJAX setup
    $.ajaxSetup({ cache: false });
    // $('[data-toggle="popover"]').popover(); 

    $("form").on('change', function(){
        $("#results").hide();
        $("#identify_alert").hide();
        // $("#identify_submit_button").show();
    });

    $("#identify_reset_button").on('click', function(){
        $("#results").hide();
        $("#identify_alert").hide();
        $('form')[0].reset();
    });

    // $("#mutations").bind("keydown", function (e) {
    //     $("#results").hide();
    //     $("#identify_alert").hide();
    // });
    function failedRequest() {
            $("#identify_alert").toggleClass('alert-warning', false);
            $("#identify_alert").toggleClass('alert-success', false);
            $("#identify_alert").toggleClass('alert-danger', true);
            $("#identify_alert_message").text("Timeout error. Try again or submit a smaller data sample.");
            $("#identify_alert").show();
            $("#identify_loading_button").hide();
            $("#identify_submit_button").show();
            return;
    }

    function updateResults(data, status){
        // console.log(data);

        if (data.hasOwnProperty('error')) {
            $("#identify_alert").toggleClass('alert-warning', false);
            $("#identify_alert").toggleClass('alert-success', false);
            $("#identify_alert").toggleClass('alert-danger', true);
            $("#identify_alert_message").text(data.error);
            $("#identify_alert").show();
            $("#identify_loading_button").hide();
            $("#identify_submit_button").show();
            return;
        } else {
            if (data.hasOwnProperty('warning')) {
                $("#identify_alert").toggleClass('alert-danger', false);
                $("#identify_alert").toggleClass('alert-success', false);
                $("#identify_alert").toggleClass('alert-warning', true);
                $("#identify_alert_message").text(data.warning);
                $("#identify_alert").show();
                // $("#identify_alert").fadeTo(1000, 1).fadeOut(600);
            }
            if (data.hasOwnProperty('success')) {
                $("#identify_alert").toggleClass('alert-danger', false);
                $("#identify_alert").toggleClass('alert-warning', false);
                $("#identify_alert").toggleClass('alert-success', true);
                $("#identify_alert_message").text(data.success);
                $("#identify_alert").show();
                // $("#identify_alert").fadeTo(1000, 1).fadeOut(600);
            }
            // FILL IN mutations input
            // $("#mutations").val(data.input);
            query_profile = data.query_formatted;
        
            drawBarplotJSON("#identify_query", data.query, 850, 300);
            $("#bt_table_type").bootstrapTable('load', data.cancer_type);
            $("#bt_table_site").bootstrapTable('load', data.primary_site);
            $("#bt_table_decomposition_A_1").bootstrapTable('load', data.decomposition_A_1);
            $("#bt_table_decomposition_B_1").bootstrapTable('load', data.decomposition_B_1);
            $("#bt_table_decomposition_cosmic_1").bootstrapTable('load', data.decomposition_cosmic_1);
            $("#bt_table_decomposition_A_2").bootstrapTable('load', data.decomposition_A_2);
            $("#bt_table_decomposition_B_2").bootstrapTable('load', data.decomposition_B_2);
            $("#bt_table_decomposition_cosmic_2").bootstrapTable('load', data.decomposition_cosmic_2);
            drawClusteringPlot("#identify_map", data.map2d, data.points);
            
            // $("#residuals_A").text(data.residuals_A);
            // $("#residuals_B").text(data.residuals_B);
            // $("#residuals_cosmic").text(data.residuals_cosmic);

            $("#bt_table_similar_samples").bootstrapTable('load', data.similar_samples);
            // showSimilarSamples("#similar_samples", data.similar_samples);
            ////////////
            $("#identify_loading_button").hide();
            $("#identify_submit_button").show();
            $("#results").show();
            scrollToAnchor("results");
        }
    }

    $("form").submit(function(e){
        e.preventDefault();
        $("#identify_submit_button").hide();
        $("#identify_loading_button").show();
        $("#identify_alert").hide();
        $("#results").hide();

        // Get the selected files from the input and store in a formData object
        var fileSelect = document.getElementById('inputFile');
        var files = fileSelect.files;
        var formData = new FormData();

        if (files.length == 1) {
            var file = files[0];
            formData.append('datafile', file, file.name);
            // console.log("Append file " +  file.name);
        }
        formData.append('mutations', $("#mutations").val());
        formData.append('format', '');  // $("#format").val());
        formData.append('mode', $("#mode").val());
        formData.append('assembly', $("#assembly").val());
        // console.log($("#assembly").val());
        formData.append('classification_method', $("#classification_method").val());

        // console.log(formData);
        // console.log(_csrf_token);

        $.ajax({
            url: url_api_identify,
            data: formData,
            processData: false,
            type: 'POST',
            contentType: false,
            // contentType: 'multipart/form-data', 
            // mimeType: 'multipart/form-data',
            beforeSend: function(xhr){
                xhr.setRequestHeader('X-CSRFToken', _csrf_token); return true;
            },
            // headers: { 'X-CSRFToken': _csrf_token,
            //             'contentType': 'multipart/form-data'},
            success: updateResults,
            failure: failedRequest
        });
    });

    ////////////////////////////////////////////////////////////////////

    var btn_png = document.getElementById('save_png');
    var btn_svg = document.getElementById('save_svg');
    var btn_data = document.getElementById('save_data');

    var svg = document.getElementById('identify_query');
    // console.log(svg);
    var canvas = document.querySelector('canvas');

    function triggerDownload (imgURI, extension) {
      var evt = new MouseEvent('click', {
        view: window,
        bubbles: false,
        cancelable: true
      });

      var a = document.createElement('a');
      a.setAttribute('download', 'profile.'+extension);
      a.setAttribute('href', imgURI);
      a.setAttribute('target', '_blank');

      a.dispatchEvent(evt);
    }


    function destroyClickedElement(event)
    {
        document.body.removeChild(event.target);
    }


    btn_svg.addEventListener('click', function () {
      var data = (new XMLSerializer()).serializeToString(svg);
      var svgBlob = new Blob([data], {type: 'image/svg+xml;charset=utf-8'});
      var downloadLink = document.createElement("a");
      downloadLink.download =  "profile.svg";
      downloadLink.innerHTML = "Download File";
      if (window.webkitURL != null)
      {
          // Chrome allows the link to be clicked
          // without actually adding it to the DOM.
          downloadLink.href = window.webkitURL.createObjectURL(svgBlob);
      }
      else
      {
          // Firefox requires the link to be added to the DOM
          // before it can be clicked.
          downloadLink.href = window.URL.createObjectURL(svgBlob);
          downloadLink.onclick = destroyClickedElement;
          downloadLink.style.display = "none";
          document.body.appendChild(downloadLink);
      }

      downloadLink.click();
    });


    btn_data.addEventListener('click', function () {
        var textToWrite = query_profile;
        var textFileAsBlob = new Blob([textToWrite], {type:'text/plain'});
        var downloadLink = document.createElement("a");
        downloadLink.download =  "profile.txt";
        downloadLink.innerHTML = "Download File";
        if (window.webkitURL != null)
        {
            // Chrome allows the link to be clicked
            // without actually adding it to the DOM.
            downloadLink.href = window.webkitURL.createObjectURL(textFileAsBlob);
        }
        else
        {
            // Firefox requires the link to be added to the DOM
            // before it can be clicked.
            downloadLink.href = window.URL.createObjectURL(textFileAsBlob);
            downloadLink.onclick = destroyClickedElement;
            downloadLink.style.display = "none";
            document.body.appendChild(downloadLink);
        }

        downloadLink.click();
    });

    btn_png.addEventListener('click', function () {        
        var svgstring = (new XMLSerializer()).serializeToString(svg);
        svgstring = window.btoa(unescape(encodeURIComponent(svgstring)));
        console.log(svgstring);
        var data = "data:image/svg+xml;base64," + svgstring;
        
        var ctx = canvas.getContext('2d');
        var img = new Image(2350, 630);
        
        img.onload = function () {
            ctx.drawImage(img, 0, 0);
            canvas.toBlobHD(function(blob) {
                saveAs(blob, "profile.png");
            });
        };
        img.src = data;
    });


});
