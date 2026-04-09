

function drawBarplot(plot_id, url_signature, max_width, max_height, onclick, show_legend, negative) {
    d3.json(url_signature, function(error, json) {
        drawBarplotJSON(plot_id, json.data, max_width, max_height, onclick, show_legend, negative);
    });
}

function drawBarplotJSON(plot_id, json, max_width, max_height, onclick, show_legend, negative, show_context) {

    d3.selectAll(plot_id + " > *").remove();

    max_width = typeof max_width !== 'undefined' ? max_width : 960;
    max_height = typeof max_height !== 'undefined' ? max_height : 400;
    show_legend = typeof show_legend !== 'undefined' ? show_legend : true;
    negative = typeof negative !== 'undefined' ? negative : false;
    show_context = typeof show_context !== 'undefined' ? show_context : true;

    var labels_offset = 0;
    if (show_context) {
        labels_offset = 23;
    }

    var margin = {top: 20, right: 75, bottom: 40 + labels_offset, left: 60},
        width = max_width - margin.left - margin.right,
        height = max_height - margin.top - margin.bottom;

    var x0 = d3.scale.ordinal()
        .rangeRoundBands([0, width], 0.07);

    var x01 = d3.scale.ordinal()
        .rangeRoundBands([0, width], 0.07);

    var x1 = d3.scale.ordinal();

    var y = d3.scale.linear()
        .range([height, 0]);
    
    // var y = d3.scale.linear()
    //     .range([height, ]);

    // var color = d3.scale.ordinal()
    //     .range(["#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00"]);
    // var color = d3.scale.category20c();
    var color = d3.scale.category10();

    var xAxis = d3.svg.axis()
        .scale(x0)
        .outerTickSize(0)
        .tickFormat(function(m){
            return m[0] + " → " + m[1];
          })
        .orient("bottom");

    var xAxis2 = d3.svg.axis()
        .scale(x01)
        .outerTickSize(0)
        .tickFormat(function(m){
            return m[0] + " → " + m[1];
          })
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left")
        .ticks(5)
        .tickFormat(d3.format(".3"));

        // .innerTickSize(6)
        // .outerTickSize(1)

    function customYAxis(g) {
      g.call(yAxis);
      g.selectAll("path").attr("fill", "none").attr("stroke-width", 1).attr("stroke", "black");
      g.selectAll(".tick line").attr("stroke-width", 0.5).attr("stroke", "black");
    }

    function customXAxis(g) {
      g.call(xAxis);
      g.select(".domain").remove();
      g.selectAll(".tick text").attr("font-size", 11);
    }

    function customXAxis2(g) {
      g.call(xAxis2);
      g.select(".domain").remove();
      g.selectAll(".tick text").attr("font-size", 11);
      //.attr("font-family", "Courier, monospace").
    }

    var tooltip = d3.select("#histogram-tooltip");

    // var svg = d3.select(plot_id)
    //     .attr("width", width + margin.left + margin.right)
    //     .attr("height", height + margin.top + margin.bottom)
    //   .append("g")
    //     .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    var svg = d3.select(plot_id)
        .attr("preserveAspectRatio", "xMaxYMax meet")
        .attr("viewBox", "0 0 " + (width + margin.left + margin.right) + " " + (height + margin.top + margin.bottom))
        .style("width", "100%")
        .style("height", "100%");

    svg.append("g")
      .append("rect")
        .attr("class", "background")
        .attr("x", 0)
        .attr("y", 0)
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .style("fill", "white");

    svg = svg.append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    $(plot_id).parent().css("height", (height + margin.top + margin.bottom) + "px");

    var data = json;

    var complementary_mutation = {};
    var complementary_context = {};
    for (i = 0, len = data.length; i < len; i++) {
      d = data[i];
      complementary_mutation[d.mutation] = d.compl_mutation;
      complementary_context[d.context] = d.compl_context;
    }

    var mutations = d3.set(data.map(function(d) { return d.mutation; })).values().sort();
    var contexts = d3.set(data.map(function(d) { return d.context; })).values().sort();
    x0.domain(mutations);
    x01.domain(mutations.map(function(d) { return complementary_mutation[d]; }));
    x1.domain(contexts).rangeRoundBands([0, x0.rangeBand()]);
    y.domain([0, d3.max(data, function(d) { return parseFloat(d.freq); }) ]);

    // initialize colors by addressing them in a predefined (sorted) order
    // var colors = contexts.map(color);
    var colors = mutations.map(color);
    // console.log(colors);

    svg.selectAll(".grid-block")
      .data(x0.range())
    .enter().append("rect")
      .attr("x", function(x){ return x; })
      .attr("width", x0.rangeBand())
      .attr("height", height + margin.bottom)
      // .style("fill",  "white");
      .style("fill", function(x, i) { if (i%2) return "rgb(96%, 96%, 96%)"; else return "rgb(96%, 100%, 100%)"; });
      // .style("fill", function(x, i) { if (i%2) return "rgb(97%, 97%, 97%)"; else return "rgb(97%, 100%, 100%)"; });
 
    svg.selectAll(".lines")
      .data(x0.range())
    .enter().append("line")
      .attr("x1", function(x){ return x; })
      .attr("y1", height + 2)
      .attr("x2", function(x){ return x + x0.rangeBand(); })
      .attr("y2", height + 2)
      .attr("stroke-width", 1.5)
      // .attr("opacity", 1)
      .attr("stroke", "black");

    svg.append("g")
        .attr("class", "x invisible-axis")
        // .attr("visibility", "hidden")
        .attr("transform", "translate(0," + (height + labels_offset + 3) + ")")
        .call(customXAxis);

    svg.append("g")
        .attr("class", "x invisible-axis")
        // .attr("visibility", "hidden")
        .attr("transform", "translate(0," + (height + labels_offset + 14) + ")")
        .call(customXAxis2);

    svg.append("g")
        .attr("transform", "translate(-15, 0)")
        .attr("class", "y axis")
        .call(customYAxis)
      .append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 5)
        .attr("dy", ".71em")
        .attr("font-size", "12px")
        .style("text-anchor", "end")
        .text("Frequency");

    // context labels for each mutation bar
    if (show_context) {
        var mutation_labels = svg.selectAll(".mutations")
            .data(mutations)
          .enter().append("g")
            .attr("class", "mutation-labels")
            .attr("transform", function(d) { return "translate(" + x0(d) + ", " + (height + 5) + ")"; });
        mutation_labels.selectAll(".mutation-labels")
            .data(function(mut_group) { return data.filter(function(z) { return z.mutation == mut_group;} ); })
            .enter().append("text")
                .attr("x", 0)
                .attr("font-family", "Courier, monospace")
                .attr("y",  function(d) { return x1(d.context); })
                .attr("font-size", 1.3 * x1.rangeBand())
                .text(function(d){ return d.context[0] + d.mutation[0] + d.context[1]; });

        mutation_labels.selectAll("text")  
            .style("text-anchor", "end")
            .attr("dx", "-1em")
            .attr("dy", "0.8em")
            .attr("transform", "rotate(-90)");      
    }

    var mutations = svg.selectAll(".mutations")
        .data(mutations)
      .enter().append("g")
        .attr("class", "mutations")
        .attr("transform", function(d) { return "translate(" + x0(d) + ",0)"; });

    // console.log(data);
    mutations.selectAll(".mutations")
        .data(function(mut_group) { return data.filter(function(z) { return z.mutation == mut_group;} ); })
        // .data(data)
      .enter().append("rect")
        .attr("class", "mutations")
        .attr("width", x1.rangeBand())
        .attr("x", function(d) { return x1(d.context); })
        .attr("y", function(d) { return y(parseFloat(d.freq)); })
        .attr("height", function(d) { return Math.max(height - y(parseFloat(d.freq)) - 1, 0); })
        .style("fill", function(d) { 
          // return color(d.context);
          return color(d.mutation);
        })
        .attr("stroke-width", 0.5)
        .attr("stroke", "black")
      .on("click", function(d) {
            // console.log("CLICK");
            if (onclick) {
              var selected = d3.select(this).classed("bar-selected");
              // remove all previously marked bars
              d3.select(".bar-selected").classed("bar-selected", false);
              // invert selector
              selected = !selected;
              d3.select(this).classed("bar-selected", selected);
              return onclick(d, selected);
            }
          })
      .on("mouseover", function(d) {
                  // console.log(d);
                  tooltip.transition()        
                      .duration(200)      
                      .style("opacity", .97);
                  tooltip.style("left", (d3.event.pageX) + "px").style("top", (d3.event.pageY - 28) + "px");
                  d3.select('text#dna1').text(d.context[0]);
                  d3.select('text#dna2').text(d.mutation[0]);
                  d3.select('text#dna3').text(d.context[1]);

                  d3.select('text#cdna1').text(d.compl_context[1]);
                  d3.select('text#cdna2').text(d.compl_mutation[0]);
                  d3.select('text#cdna3').text(d.compl_context[0]);

                  d3.select('text#mut').text(d.mutation[1]);
                  d3.select('text#cmut').text(d.compl_mutation[1]);

                  d3.select('#ncases').html(d.count);  
                  })                  
      .on("mouseout", function(d) {       
                  tooltip.transition()        
                      .duration(400)      
                      .style("opacity", 0);  
              });

    if (false && show_legend) {
          var legend = svg.selectAll(".legend")
              .data(contexts.slice())
            .enter().append("g")
              .attr("class", "legend")
              .attr("transform", function(d, i) { return "translate(0," + (i * 16 + 10) + ")"; });

          legend.append("rect")
              .attr("x", width+margin.right - 30)
              .attr("width", 12)
              .attr("height", 12)
              .style("fill", color);

          legend.append("text")
              .attr("x", width+margin.right - 35)
              .attr("y", 6)
              .attr("dy", 4)
              .attr("font-family", "Courier")
              .attr("font-size", "12px")
              .style("text-anchor", "end")
              .text(function(d) { return d[0] + '⋅x⋅' + d[1]; });

          svg.append("text")
              .attr("x", width - 5)
              .attr("y", 0)
              .attr("font-size", "10px")
              .text("5'⋅x⋅3' context");
    }

} // drawBarplot




function drawMutabilityPlot(handle, url) {
    d3.json(url, function(error, json) {
        drawMutabilityPlotJSON(handle, json);
    });
}


function drawMutabilityPlotJSON(handle, json) {
    
    max_width = typeof max_width !== 'undefined' ? max_width : 720;
    max_height = typeof max_height !== 'undefined' ? max_height : 500;
    show_legend = typeof show_legend !== 'undefined' ? show_legend : true;
    negative = typeof negative !== 'undefined' ? negative : false;

    var margin = {top: 5, right: 5, bottom: 40, left: 40},
        width = max_width - margin.left - margin.right,
        height = max_height - margin.top - margin.bottom;

    var y = d3.scale.ordinal()
        .rangeRoundBands([0, height], 0.07);

    // var x01 = d3.scale.ordinal()
    //     .rangeRoundBands([0, width], 0.07);

    // var x1 = d3.scale.ordinal();

    var x = d3.scale.linear()
        .range([0, width]);
    
    // var y = d3.scale.linear()
    //     .range([height, ]);

    // var color = d3.scale.ordinal()
    //     .range(["#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00"]);
    // var color = d3.scale.category10();
    // 
    //  var color = d3.scale.category10();
    var color = function(n){
        // #ad494a #8ca252
        var colors = ["#5ab4ac", "#d6616b", ];
        // var colors = ["#d0743c", "#7b6888"];
        return d3.scale.ordinal().range(colors);
    }();

    // var xAxis = d3.svg.axis()
    //     .scale(x0)
    //     .tickFormat(function(m){
    //         return m[0] + " → " + m[1];
    //       })
    //     .orient("bottom");

    // var xAxis2 = d3.svg.axis()
    //     .scale(x01)
    //     .tickFormat(function(m){
    //         return m[0] + " → " + m[1];
    //       })
    //     .orient("bottom");

    // var yAxis = d3.svg.axis()
    //     .scale(y)
    //     .orient("left")
    //     .tickFormat(d3.format(".3"));

    // var tooltip = d3.select("#histogram-tooltip");

    // var svg = d3.select(plot_id)
    //     .attr("width", width + margin.left + margin.right)
    //     .attr("height", height + margin.top + margin.bottom)
    //   .append("g")
    //     .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    var svg = d3.select(handle)
        .attr("preserveAspectRatio", "xMaxYMax meet")
        .attr("viewBox", "0 0 " + max_width + " " + max_height)
        .style("width", max_width)
        .style("height", max_height);
    $(handle).parent().css("height", (max_height) + "px");
      // .append("g")
      //   .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

      var data = json.data;
      // console.log(data);
      // console.log(svg);

      // var complementary_mutation = {};
      // var complementary_context = {};
      // for (i = 0, len = data.length; i < len; i++) {
      //   d = data[i];
      //   complementary_mutation[d.mutation] = d.compl_mutation;
      //   complementary_context[d.context] = d.compl_context;
      // }

      var trinucleotides = d3.set(data.map(function(d) { return d.trinucleotide; })).values().sort();
      // var contexts = d3.set(data.map(function(d) { return d.context; })).values().sort();
      y.domain(trinucleotides);
      // x01.domain(mutations.map(function(d) { return complementary_mutation[d]; }));
      // x1.domain(contexts).rangeRoundBands([0, x0.rangeBand()]);
      x.domain([0, d3.max(data, function(d) { return parseFloat(d.mutability); }) ]);

      // initialize colors by addressing them in a predefined (sorted) order
      // var colors = contexts.map(color);
      // var colors = mutations.map(color);
      // console.log(colors);
    
      // svg.selectAll(".grid-block")
      //   .data(x0.range())
      // .enter().append("rect")
      //   .attr("x", function(x){ return x; })
      //   .attr("width", x0.rangeBand())
      //   .attr("height", height + margin.bottom)
      //   // .style("fill",  "white");
      //   .style("fill", function(x, i) { if (i%2) return "rgb(96%, 96%, 96%)"; else return "rgb(96%, 100%, 100%)"; });
      //   // .style("fill", function(x, i) { if (i%2) return "rgb(97%, 97%, 97%)"; else return "rgb(97%, 100%, 100%)"; });
   
      // svg.selectAll(".lines")
      //   .data(x0.range())
      // .enter().append("rect")
      //   .attr("x", function(x){ return x; })
      //   .attr("y", height + 1)
      //   .attr("width", x0.rangeBand())
      //   .attr("height", 2)
      //   .style("fill",  "black");

      // svg.append("g")
      //     .attr("class", "x invisible-axis")
      //     .attr("transform", "translate(0," + (height + 7) + ")")
      //     .call(xAxis.tickSize(0, 0, 0));

      // svg.append("g")
      //     .attr("class", "x invisible-axis")
      //     .attr("transform", "translate(0," + (height + 17) + ")")
      //     .call(xAxis2);

      // svg.append("g")
      //     .attr("transform", "translate(-15, 0)")
      //     .attr("class", "y axis")
      //     .call(yAxis)
      //   .append("text")
      //     .attr("transform", "rotate(-90)")
      //     .attr("y", 6)
      //     .attr("dy", ".71em")
      //     .style("text-anchor", "end");
      //     // .text("Number of Mutations");

      // var mutations = svg.selectAll(".mutations")
      //     .data(mutations)
      //   .enter().append("g")
      //     .attr("class", "mutations")
      //     .attr("transform", function(d) { return "translate(" + x0(d) + ",0)"; });

      // console.log(data);
      var trinucl = svg.selectAll(".trinucleotide")
          // .data(function(mut_group) { return data.filter(function(z) { return z.mutation == mut_group;} ); })
          .data(data)
        .enter().append("rect")
          .attr("class", "trinucleotide")
          .attr("height", y.rangeBand())
          .attr("y", function(d) { return y(d.trinucleotide); })
          .attr("x", function(d) { return margin.left; })
          .attr("width", function(d) { return Math.max(x(parseFloat(d.mutability)) - 1, 0); })
          .style("fill", function(d) { 
            // return color(d.context);
              return color(d.trinucleotide.substring(1, 2));
          })
          // .style("fill", function(d) {
          //    return "#f00";
          // })
          .attr("stroke-width", 0.5)
          .attr("stroke", "black")
        // .on("click", function(d) {
        //       // console.log("CLICK");
        //       if (onclick) {
        //         var selected = d3.select(this).classed("bar-selected");
        //         // remove all previously marked bars
        //         d3.select(".bar-selected").classed("bar-selected", false);
        //         // invert selector
        //         selected = !selected;
        //         d3.select(this).classed("bar-selected", selected);
        //         return onclick(d, selected);
        //       }
        //     })
        // .on("mouseover", function(d) {
        //             // console.log(d);
        //             tooltip.transition()        
        //                 .duration(200)      
        //                 .style("opacity", .97);
        //             tooltip.style("left", (d3.event.pageX) + "px").style("top", (d3.event.pageY - 28) + "px");
        //             d3.select('text#dna1').text(d.context[0]);
        //             d3.select('text#dna2').text(d.mutation[0]);
        //             d3.select('text#dna3').text(d.context[1]);

        //             d3.select('text#cdna1').text(d.compl_context[1]);
        //             d3.select('text#cdna2').text(d.compl_mutation[0]);
        //             d3.select('text#cdna3').text(d.compl_context[0]);

        //             d3.select('text#mut').text(d.mutation[1]);
        //             d3.select('text#cmut').text(d.compl_mutation[1]);

        //             d3.select('#ncases').html(d.count);  
        //             })                  
        // .on("mouseout", function(d) {       
        //             tooltip.transition()        
        //                 .duration(400)      
        //                 .style("opacity", 0);  
        //         });

        var legend = svg.selectAll(".legend")
            .data(data)
          .enter().append("g")
            .attr("class", "legend")
            .attr("transform", function(d, i) { return "translate("+ margin.left +"," + y(d.trinucleotide) + ")"; });

        legend.append("text")
            .attr("x", -3)
            .attr("y", 0)
            .attr("dy", 7)
            .attr("font-family", "Courier")
            .attr("font-size", "10px")
            .style("text-anchor", "end")
            .text(function(d) { return d.trinucleotide; });


        var xAxis = d3.svg.axis()
            .scale(x)
            .orient("bottom")
            //.tickFormat(d3.format("e"))
            .ticks(5);
        svg.append("g")
            .attr("class", "x axis")
            .attr("transform", "translate(" + margin.left + "," + (height) + ")")
            .call(xAxis)
        .append("text")
          // .attr("transform", "rotate(-90)")
          .attr("x", width - 150)
          .attr("dy", "40px")
          .attr("font-size", "12px")
          .style("text-anchor", "end")
          .text("DNA mutability (expected number of mutations per exome Megabase)");

}


function show_signature_modal(id, label, header) {
    if (typeof header !== 'undefined' && header.length > 0) {
        $("#signature_modal_header").text('Mutational profile/signature for ' + header);
    } else {
        $("#signature_modal_header").text("");
    }

    d3.selectAll("#hist_modal > *").remove();
    drawBarplot("#hist_modal", encodeURI(url_signature_freqs.replace('0', id)), max_width=750, max_height=280, click='undefined', show_legend=true);
    $("#hist_modal_loading").text(label);
    $("#hist_modal_link").html('<a href="' + encodeURI(url_signature.replace('0', id)) + '">Click here to explore the signature and the corresponding mutations</a>');
    $("#signatureModal").modal('show');
} // show_signature_modal



function show_mutability_modal(id, label, header) {
    if (typeof header !== 'undefined' && header.length > 0) {
        $("#signature_modal_header").text('Mutability profile for ' + header);
    } else {
        $("#signature_modal_header").text("");
    }

    d3.selectAll("#hist_modal > *").remove();
    // drawBarplot("#hist_modal", encodeURI(url_signature_freqs.replace('0', id)), max_width=750, max_height=280, click=undefined, show_legend=true);
    drawMutabilityPlot("#hist_modal", encodeURI(url_mutability.replace('0', id)), max_width=750, max_height=280, click='undefined', show_legend=true)
    $("#hist_modal_loading").text(label);
    $("#hist_modal_link").html('<a href="' + encodeURI(url_signature.replace('0', id)) + '">Click here to explore the corresponding mutational profile/signature</a>');
    $("#signatureModal").modal('show');
} // show_signature_modal


function drawExposureMatrix(handle, data) {
      var margin = {top: 0, right: 250, bottom: 100, left: 20},
          width = 1100 - margin.left - margin.right,
          height = 500 - margin.top - margin.bottom;

      // console.log(data);

      var x = d3.scale.ordinal()
          .rangeRoundBands([0, width]);

      var y = d3.scale.linear()
          .rangeRound([height, 0]);

      // var z = d3.scale.category20();
      var z = d3.scale.ordinal()
        .range(["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99"]);

      var xAxis = d3.svg.axis()
          .scale(x)
          .orient("bottom");

      var yAxis = d3.svg.axis()
          .scale(y)
          .orient("right");

      var svg = d3.select(handle)
          .attr("width", width + margin.left + margin.right)
          .attr("height", height + margin.top + margin.bottom)
        .append("g")
          .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

      var signatures = data.columns;
      // console.log(signatures);
      var layers = d3.layout.stack()(signatures.map(function(c, i) {
          return data.rows.map(function(d) {
              return {x: d.sample, y: d[String(i)]};
          });
      }));

      x.domain(layers[0].map(function(d) { return d.x; }));
      y.domain([0, d3.max(layers[layers.length - 1], function(d) { return d.y0 + d.y; })]).nice();

      var layer = svg.selectAll(".layer")
          .data(layers)
        .enter().append("g")
          .attr("class", "layer")
          .style("fill", function(d, i) { return z(i); })

      // .append("svg:title").text(function(d, i) { console.log(i); });

      layer.selectAll("rect")
          .data(function(d) { return d; })
        .enter().append("rect")
          .attr("x", function(d) { return x(d.x); })
          .attr("y", function(d) { return y(d.y + d.y0); })
          .attr("height", function(d) { return y(d.y0) - y(d.y + d.y0); })
          .attr("width", x.rangeBand() );

      svg.append("g")
          .attr("class", "axis axis--x")
          .attr("transform", "translate(0," + height + ")")
          .call(xAxis)
          .selectAll("text")  
            .style("text-anchor", "end")
            .attr("dx", "-.8em")
            .attr("dy", ".15em")
            .attr("transform", function(d) {
                return "rotate(-65)" 
                });
      // svg.append("g")
      //     .attr("class", "axis axis--y")
      //     .attr("transform", "translate(" + width + ",0)")
      //     .call(yAxis);

      var legend = svg.append("g")
          .style("font-family", "sans-serif")
          .style("font-size", 10)
          .style("text-anchor", "end")
        .selectAll("g")
        .data(signatures.reverse())
        .enter().append("g")
          .attr("transform", function(d, i) { return "translate(0," + i * 20 + ")"; });

      legend.append("rect")
          .attr("x", width + margin.right - 19)
          .attr("width", 19)
          .attr("height", 19)
          .style("fill", function(d, i) {
            // console.log(d, i, signatures.length - i - 1);
            return z(signatures.length - i - 1);
          });

      legend.append("text")
          .attr("x", width + margin.right - 24)
          .attr("y", 9.5)
          .attr("dy", "0.32em")
          .text(function(d, j) {
              if (signatures.length > 10) {
                  var s = 0.0;

                  for (var i = 0; i < data.rows.length; i++) {
                      s += data.rows[i][j+1];
                  }

                  return d.split(" ").slice(0, 2).join(" ") + " [" + s.toFixed(2).toString() + "]";
              } else {
                  return d;
              }
          }).append("svg:title").text(function(d) { return d; });
}


function drawSamplesMatrix(svg_id, json){
    var spacer = 0;
    var block_spacer = 22;
    var XOffset = 80;
    var YOffset = 50;

    // var colorScale = d3.scale.linear().domain([0, 0.01, 0.1]).range(["white", "orange", "red"]);    // yellow - red
    var opacity = d3.scale.sqrt().domain([0.01, 0.03]).range([0.4, 1.0]).clamp(true);

    var data = json.data.matrix;
    var samples = json.data.rows;
    var Y = samples.length;
    var X = json.data.cols.length;

    // console.log(json.data);

    var w = 7; //(960.0 - 60 - 75 - block_spacer * 5) / 96.0 // 6.8;
    var h = 3;
    var size = d3.scale.linear().domain([10, 500]).rangeRound([5, 1]).clamp(true);
    h = size(Y); // scale down with the number of samples
    // h = w;

    //var mutations = d3.set(data.map(function(d) { return d.mutation; })).values().sort();
    var mutations = json.data.mutations;
    var compl_mutations = json.data.compl_mutations;
    var contexts = json.data.contexts;
    var compl_contexts = json.data.compl_contexts;

    // console.log(contexts)
    
    var color = d3.scale.category10();
    var colors = mutations.map(color);
    // console.log(mutations);

    // console.log(X, Y);
    var width = 960 - XOffset; // w*X + spacer*(X-1) + block_spacer * mutations.length;
    var height = h*Y + spacer*(Y-1);
    // console.log(svg_id);

    // max_width = typeof max_width !== 'undefined' ? max_width : 960;
    // max_height = typeof max_height !== 'undefined' ? max_height : 400;

    var svg = d3.select(svg_id)
        .attr("preserveAspectRatio", "xMaxYMax meet")
        .attr("viewBox", "0 0 " + (width + XOffset) + " " + (height + YOffset))
        .style("width", "100%")
        .style("height", "100%");

    $(svg_id).parent().css("height", (height + YOffset) + "px");
      // .append("g")
      //   .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    // var svg = d3.select(svg_id)
    //     .attr("preserveAspectRatio", "xMaxYMax meet")
    //     .attr("viewBox", "0 0 " + (width  + XOffset) + " " + (height + YOffset))
    //     .attr("width", max_width)
    //     .attr("height", height + YOffset);
    //     // .attr("width", "100%")
    //     // .attr("height", "100%");
    var tooltip = d3.select("#helper-tooltip");

    var x0 = d3.scale.ordinal()
        .domain(mutations)
        .rangeRoundBands([XOffset, width], 0.07);

    var x0c = d3.scale.ordinal()
        .domain(compl_mutations)
        .rangeRoundBands([XOffset, width], 0.07);

    var x1 = d3.scale.ordinal()
        .domain(contexts)
        .rangeRoundBands([0, x0.rangeBand()]);

    //////////////////////////////////
    svg.append("g")
      .append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 6)
        .attr("dy", ".71em")
        .attr("font-size", "12px")
        .style("text-anchor", "end")
        .text("Cancer samples");

    // var mutations_blocks = svg.selectAll(".mutations")
    //     .data(mutations)
    //   .enter().append("g")
    //     .attr("class", "mutations")
    //     .attr("transform", function(d) { return "translate(" + x0(d) + ",0)"; });
    //////////////////////////////////
    // if (show_context) {
    //     var mutation_labels = svg.selectAll(".mutations")
    //         .data(mutations)
    //       .enter().append("g")
    //         .attr("class", "mutation-labels")
    //         .attr("transform", function(d) { return "translate(" + x0(d) + ", " + (height + 0) + ")"; });
    //     mutation_labels.selectAll(".mutation-labels")
    //         .data(function(mut_group) { return data.filter(function(z) { return z.mutation == mut_group;} ); })
    //         .enter().append("text")
    //             .attr("x", 0)
    //             .attr("font-family", "Courier, monospace")
    //             .attr("y",  function(d) { return x1(d.context); })
    //             .attr("font-size", 1.3 * x1.rangeBand())
    //             .text(function(d){ return d.context[0] + d.mutation[0] + d.context[1]; });

    //     mutation_labels.selectAll("text")  
    //         .style("text-anchor", "end")
    //         .attr("dx", "-1em")
    //         .attr("dy", "0.8em")
    //         .attr("transform", "rotate(-90)");      
    // }







    // var mutation_labels = svg.selectAll(".mutations")
    //     .data(contexts)
    //   .enter().append("g")
    //     .attr("class", "mutation-labels")
    //     .attr("transform", function(d, i) { console.log(d,  x1(d)); return "translate(" + x1(d) + ", " + (100) + ")"; });

    // mutation_labels.selectAll(".mutation-labels")
    //     .data(mutations)
    //     .enter().append("text")
    //         .attr("x", 0)
    //         .attr("font-family", "Courier, monospace")
    //         .attr("y",  function(d) { return x1(d.context); })
    //         .attr("font-size", 1.3 * x1.rangeBand())
    //         .text(function(d, i){
    //               return d[0] + mutations[i][0] + d[1]; 
    //         });
    // mutation_labels.selectAll("text")  
    //     .style("text-anchor", "end")
    //     .attr("dx", "-1em")
    //     .attr("dy", "0.8em")
    //     .attr("transform", "rotate(-90)");      


    //////////////////////////////////

    svg.append("g")
        .attr("class", "x invisible-axis")
        .attr("transform", "translate(0," + 7 + ")")
        .call(
              d3.svg.axis()
              .scale(x0)
              .tickFormat(function(m){
                  return m[0] + " → " + m[1];
                })
              .orient("top").tickSize(0, 0, 0)
        );

    svg.append("g")
        .attr("class", "x invisible-axis")
        .attr("transform", "translate(0," + 27 + ")")
        .call(
              d3.svg.axis()
              .scale(x0c)
              .tickFormat(function(m){
                  return m[0] + " → " + m[1];
                })
              .orient("top")
          );
    //////////////////////////////////


    var cluster_colors = ["#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"];
    var romans = ["I", "II", "III", "IV", "V", "VI"];
    ////////////////////////////////////
    // console.log(data);
    var plotCell = svg.selectAll(".sample-mutations")
         .data(data).enter().append("rect")
             // .attr('r', 1.1*w/2)
             // .attr('r', 1.1*size(Y)/2)
             .attr('width', w)
             .attr('height', h)
             .attr('x', function(d) {
                // console.log(d.x);
                // console.log(mutations.indexOf(d.mutation))
                return +d.x * (w + spacer) + (block_spacer * mutations.indexOf(d.mutation) / 16) + XOffset;
             })
             .attr('y', function(d) {
                return +d.y * (h + spacer) + YOffset;
             })
             .style('fill-opacity', function(d) {
                var value = +d.value || 0;
                return opacity(value);
             })
             .style('fill', function(d) {
                // var value = +d.value || 0;
                return color(d.mutation);
             });


    var clusterRow = svg.selectAll(".sample-row2")
        .data(json.data.rows).enter().append("rect")
            .attr('width', 60)
            .attr('height', h)
            .attr('x', function(d, i){
                return 960-80;
            })
            .attr('y', function(d, i){
                return i * (h+spacer) + YOffset;
            })
            // .style('fill-opacity', 0.0)
            .style('fill', function(d, i){
                var s = samples[i];
                if (typeof s.cluster !== 'undefined' && s.cluster >= 0) {
                    return cluster_colors[s.cluster];
                } else {
                    return "white";
                }
            });


    if (typeof json.data.flags !== 'undefined') {
        var flagRow = svg.selectAll(".sample-row3")
            .data(json.data.flags).enter().append("rect")
                .attr('width', 10)
                .attr('height', h)
                .attr('x', function(d, i){
                    return 960-10;
                })
                .attr('y', function(d, i){
                    return i * (h+spacer) + YOffset;
                })
                // .style('fill-opacity', 0.0)
                .style('fill', function(d, i){
                    // var f = json.data.flags[i];
                    if (typeof d !== 'undefined' && d > 0) {
                        return "black";
                    } else {
                        return "white";
                    }
                });      
    }

    var plotRow = svg.selectAll(".sample-row")
        .data(json.data.rows).enter().append("rect")
            .attr('width', width + XOffset)
            .attr('height', h)
            .attr('x', function(d, i){
                return 0;
            })
            .attr('y', function(d, i){
                return i * (h+spacer) + YOffset;
            })
            .style('fill-opacity', 0.0)
            .style('fill', '#888')
        .on("mouseover", function(d, i) {
                    var s = samples[i];
                    
                    d3.select(this).style('fill-opacity', 0.5);
                    tooltip.style("left", (d3.event.pageX) + "px").style("top", (d3.event.pageY - 28) + "px");
                    tooltip.transition()
                        .duration(200)
                        .style("opacity", .97);
                    // var html = "<b>Sample " + s.sample_name + "</b><br>";
                    var html = "Sample with ";
                    if (s.n_mutations) {
                         html += s.n_mutations + " mutations<br><br>";
                    }
                    if (s.age) {
                         html += "Age: "+ s.age + "<br>";
                    }
                    if (s.gender) {
                        if (s.gender == 'm') {
                          html += "Male <br>";
                        }
                        else if (s.gender == "f") {
                          html += "Female <br>";
                          
                        }
                    }
                    if (s.stage) {
                         html += "Stage: "+ s.stage + "<br>";
                    }
                    if (s.ethnicity) {
                         html += "Ethnicity: "+ s.ethnicity + "<br>";
                    }
                    if (s.therapy) {
                         html += "Therapy: "+ s.therapy + "<br>";
                    }
                    if (s.environment) {
                         html += "Environm.: "+ s.environment + "<br>";
                    }
                    if (s.msi) {
                         html += "MSI: "+ s.msi + "<br>";
                    }
                    if (typeof s.cluster !== 'undefined' && s.cluster >= 0)  {
                         html += "<br>Cluster: "+ romans[s.cluster] + "<br>";
                    }
                    tooltip.html(html);

                    // console.log(d);
                    // tooltip.transition()        
                    //     .duration(200)      
                    //     .style("opacity", .97);
                    // tooltip.style("left", (d3.event.pageX) + "px").style("top", (d3.event.pageY - 28) + "px");
                    // d3.select('text#dna1').text(d.context[0]);
                    // d3.select('text#dna2').text(d.mutation[0]);
                    // d3.select('text#dna3').text(d.context[1]);

                    // d3.select('text#cdna1').text(d.compl_context[1]);
                    // d3.select('text#cdna2').text(d.compl_mutation[0]);
                    // d3.select('text#cdna3').text(d.compl_context[0]);

                    // d3.select('text#mut').text(d.mutation[1]);
                    // d3.select('text#cmut').text(d.compl_mutation[1]);

                    // d3.select('#ncases').html(d.count);

              })                  
        .on("mouseout", function(d) {
                  d3.select(this).style('fill-opacity', 0.0);
                    // console.log("out");
                    // d3.select(".sample-row").style('fill-opacity', 0.0);  
                    tooltip.transition()        
                        .duration(400)      
                        .style("opacity", 0);  
                });

} //drawSamplesMatrix


