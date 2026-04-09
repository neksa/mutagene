

function draw_scatterplot(ref, json, logscale) {
    logscale = typeof logscale !== 'undefined' ? logscale : false;

    $("#scatterplot_loading").text("");
    d3.selectAll(ref + " > *").remove();
    var data = json.data;

    var dot_size = 5;
    var width = 500;
    var height = 430;
    var margin = 80;
    var svg = d3.select(ref)
        .attr("width", width + 2 * margin)
        .attr("height", height + margin)
      .append("g")
        .attr("transform", "translate(" + (margin) + "," + 10 + ")");

    // d3.select("#hist2")
    //   .attr("width", 1)
    //   .attr("height", 1);

    var xValue = function(d) { return d.v1;}, // data -> value
        xScale = d3.scale.linear().range([0, width]), // value -> display
        xMap = function(d) { return xScale(xValue(d));}, // data -> display
        xAxis = d3.svg.axis().scale(xScale).orient("bottom");

    var yValue = function(d) { return d.v2;}, // data -> value
        yScale = d3.scale.linear().range([height, 0]), // value -> display
        yMap = function(d) { return yScale(yValue(d));}, // data -> display
        yAxis = d3.svg.axis().scale(yScale).orient("left");

    var contexts = d3.set(data.map(function(d) { return d.attrib.context; })).values().sort();
    var mutations = d3.set(data.map(function(d) { return d.attrib.mut; })).values().sort();

    var color = d3.scale.category10();
    var colors = mutations.map(color);

    var min_val = Math.min(d3.min(data, yValue), d3.min(data, xValue)) - 0.01;
    var max_val = Math.max(d3.max(data, yValue), d3.max(data, xValue)) + 0.01;
    xScale.domain([min_val, max_val]);
    yScale.domain([min_val, max_val]);
    // xScale.domain([d3.min(data, xValue) - 0.0001, d3.max(data, xValue)+ 0.01] );
    // yScale.domain([d3.min(data, yValue) - 0.0001, d3.max(data, yValue)+ 0.01] );

    var log = "";
    if (logscale) {
      log = " (log scale)";
    }

    // x-axis
    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + height + ")")
        .call(xAxis)
      .append("text")
        .attr("class", "scatterplot-label")
        .attr("x", width)
        .attr("y", 35)
        .style("text-anchor", "end")
        .text(function(){
          return "Mutation frequencies in " + json.xlabel + log;
        });

    // y-axis
    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis)
      .append("text")
        .attr("class", "scatterplot-label")
        .attr("transform", "rotate(-90)")
        .attr("y", -55)
        .attr("dy", ".71em")
        // .attr("dy", "1em")
        .style("text-anchor", "end")
        .text(function(){
          return "Mutation frequencies in " + json.ylabel + log;
        });

    var tooltip = d3.select(".scatterplot-tooltip");

    svg.selectAll(".dot")
         .data(data)
           .enter().append("circle")
             .attr('class', 'scatter-dot')
             .attr('r', dot_size)
             .attr('cx', xMap)
             .attr('cy', yMap)
             .style('fill', function(d) {
                return color(d.attrib.mut);
             })
             .on("mouseover", function(d) {
                 tooltip.transition()
                      .duration(300)
                      .style("opacity", .95);
                 tooltip.html(d.attrib.mutation + "<br/> (" + xValue(d).toFixed(2)
                + ", " + yValue(d).toFixed(2) + ")")
                      .style("left", (d3.event.pageX + 5) + "px")
                      .style("top", (d3.event.pageY - 28) + "px");
             })
             .on("mouseout", function(d) {
                 tooltip.transition()
                      .duration(500)
                      .style("opacity", 0);
             });

    var show_legend = true;
    if (show_legend) {
    
      var legend = svg.selectAll(".legend")
          .data(mutations.slice())
        .enter().append("g")
          .attr("class", "legend")
          .attr("transform", function(d, i) { return "translate(0," + (i * 16 + 10) + ")"; });

      legend.append("circle")
          .attr("cx", width+margin - 10)
          .attr('r', 6)
          .attr('class', 'scatter-dot')
          .style("fill", color);

      legend.append("text")
          .attr("x", width+margin - 20)
          .attr("y", 0)
          .attr("dy", 4)
          .attr("font-family", "Courier")
          .attr("font-size", "12px")
          .style("text-anchor", "end")
          .text(function(d) { return d[0] + ' → ' + d[1]; });

      svg.append("text")
          .attr("x", width + margin - 0)
          .attr("y", 0)
          .style("text-anchor", "end")
          .text("Mutation:");
    }
}

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

function drawLogPlot(plot_id, url_signature) {

    max_width = typeof max_width !== 'undefined' ? max_width : 720;
    max_height = typeof max_height !== 'undefined' ? max_height : 280;
    show_legend = typeof show_legend !== 'undefined' ? show_legend : true;
    negative = typeof negative !== 'undefined' ? negative : false;

    var margin = {top: 20, right: 75, bottom: 50, left: 60},
        width = max_width - margin.left - margin.right,
        height = max_height - margin.top - margin.bottom;

    var x0 = d3.scale.ordinal()
        .rangeRoundBands([0, width], 0.07);

    var x01 = d3.scale.ordinal()
        .rangeRoundBands([0, width], 0.07);

    var x1 = d3.scale.ordinal();

    var y = d3.scale.linear()
        .range([height, 0]);
    
    // var color = d3.scale.ordinal()
    //     .range(["#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00"]);
    var color = d3.scale.category10();
    // var color = d3.scale.category20c();

    var xAxis = d3.svg.axis()
        .scale(x0)
        .tickFormat(function(m){
            return m[0] + " → " + m[1];
          })
        .orient("bottom");

    var xAxis2 = d3.svg.axis()
        .scale(x01)
        .tickFormat(function(m){
            return m[0] + " → " + m[1];
          })
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left")
        .tickFormat(d3.format(".3"));

    var tooltip = d3.select("#histogram-tooltip");

    var svg = d3.select(plot_id)
        .attr("preserveAspectRatio", "xMaxYMax meet")
        .attr("viewBox", "0 0 " + (width + margin.left + margin.right) + " " + (height + margin.top + margin.bottom))
        .style("width", "100%")
        .style("height", "100%")
      .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
    $(plot_id).parent().css("height",  (height + margin.top + margin.bottom) + "px");

    d3.json(url_signature, function(error, json) {
      if (error) throw error;
      var data = json.data;
      // console.log(json);

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
      y.domain([
            d3.min(data, function(d) { return parseFloat(d.freq) - 0.1 ; }),
            d3.max(data, function(d) { return parseFloat(d.freq) + 0.1 ; })
          ]);

      // initialize colors by addressing them in a predefined (sorted) order
      var colors = mutations.map(color);
    
      svg.selectAll(".grid-block")
        .data(x0.range())
      .enter().append("rect")
        .attr("x", function(x){ return x; })
        .attr("width", x0.rangeBand())
        .attr("height", height + margin.bottom)
        // .style("fill",  "white");
        .style("fill", function(x, i) { if (i%2) return "rgb(96%, 96%, 96%)"; else return "rgb(96%, 100%, 100%)"; });
        // .style("fill", function(x, i) { if (i%2) return "rgb(97%, 97%, 97%)"; else return "rgb(97%, 100%, 100%)"; });
   
      svg.append("rect")
        .attr("x", x0(0))
        .attr("y", y(0))
        .attr("width", width)
        .attr("height", 1)
        .style("fill", "black");

      svg.selectAll(".lines")
        .data(x0.range())
      .enter().append("rect")
        .attr("x", function(x){ return x; })
        .attr("y", height + 1)
        .attr("width", x0.rangeBand())
        .attr("height", 2)
        .style("fill",  "black");

      svg.append("g")
          .attr("class", "x invisible-axis")
          .attr("transform", "translate(0," + (height + 7) + ")")
          .call(xAxis.tickSize(0, 0, 0));

      svg.append("g")
          .attr("class", "x invisible-axis")
          .attr("transform", "translate(0," + (height + 17) + ")")
          .call(xAxis2);

      svg.append("g")
          .attr("transform", "translate(-15, 0)")
          .attr("class", "y axis")
          .call(yAxis)
        .append("text")
          .attr("transform", "rotate(-90)")
          .attr("y", 5)
          .attr("dy", ".71em")
          .attr("font-size", "12px")
          .style("text-anchor", "end")
          .text("Log-ratio of mutational frequencies");

      var mutations = svg.selectAll(".mutation")
          .data(data)
        .enter().append("g")
          .attr("class", "g")
          .attr("transform", function(d) { return "translate(" + x0(d.mutation) + ",0)"; });

      mutations.selectAll("rect")
          .data(function(d) { return data.filter(function(z) { return z.mutation == d.mutation;} ); })
        .enter().append("rect")
          .attr("class", "mutation-bar")
          .attr("width", x1.rangeBand())
          .attr("x", function(d) { return x1(d.context); })
          .attr("y", function(d) {
            // console.log(d.context, d.mutation, d.freq, y(parseFloat(d.freq)));
            if (d.freq > 0) {
              return y(parseFloat(d.freq)) - 1;
            } else {
              return y(0) + 2;
            }
          })
          .attr("height", function(d) {
            if (d.freq > 0) {
                return y(0) - y(parseFloat(d.freq));
            } else if (d.freq < 0) {
                return y(parseFloat(d.freq)) - y(0);
            } else {
              return 0;
            }
          })
          .style("fill", function(d) { return color(d.mutation); })
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
                        .style("opacity", .95);
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

      if (show_legend) {
      
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


    });
}



// //////////////////////////////////////////
// // MATRIX                   
// //////////////////////////////////////////
// function show_comparison_for_two(json, placeholder_id) {
//     placeholder_id = typeof placeholder_id !== 'undefined' ? placeholder_id : "#matrix_placeholder";

//     $(placeholder_id).html("");
//     console.log(json);
    
//     var label1 = json.data.columns[0];
//     var label2 = json.data.rows[0];
//     var d = json.data.matrix[0];
//     var accession1 = json.data.column_accessions[0];
//     var accession2 = json.data.row_accessions[0];

//     d3.selectAll("#plot1 > *").remove();
//     d3.selectAll("#plot2 > *").remove();
//     $("#signaturesModal").modal('show');
//     var url_signature1 = encodeURI(url_signature_freqs.replace('0', accession1));
//     var url_signature2 = encodeURI(url_signature_freqs.replace('0', accession2));
//     // Show signature 1 and signature 2 barplots
//     $("#signature1").text(label1);
//     $("#signature2").text(label2);
//     var value = +d.value || 0;
//     var method = $("#comp_for_2").val();
//     if (method == "pearson") {
//       $("#comparison_result").text("Pearson correlation coefficient = " + value + " with p-value " + d.pvalue);
//     } else if (method == "cos") {
//       $("#comparison_result").text("Cosine coefficient = " + value);
//     } else if (method == "wilcox") {
//       $("#comparison_result").text("Wilcoxon correlation coefficient = " + value + " with p-value " + d.pvalue);
//     } else if (method == "euclidean") {
//       $("#comparison_result").text("Euclidean distance = " + value);
//     }
//     // function capitalizeFirstLetter(string) {
//     //     return string.charAt(0).toUpperCase() + string.slice(1);
//     // }
//     drawBarplot("#plot1", url_signature1, max_width=580, max_height=300, click=undefined, show_legend=true);
//     drawBarplot("#plot2", url_signature2, max_width=580, max_height=300, click=undefined, show_legend=false);
// }


function draw_matrix(json, placeholder_id, method) {
    // console.log('OK');
    placeholder_id = typeof placeholder_id !== 'undefined' ? placeholder_id : "#matrix_placeholder";
    method = typeof method !== 'undefined' ? method: $("#select3").val();

    var rows = json.data.rows;
    var cols = json.data.columns;
    var data = json.data.matrix;
    var col_accessions = json.data.column_accessions;
    var row_accessions = json.data.row_accessions;
    var mutation_type_labels = json.data.mutation_types;
    //height and width of each row in the heatmap
    //var h = 17; // 15
    var h = 13; // 15
    var w = h;
    // var w = 12.9;
    var spacer = 1;

    var x_labels_max = d3.max(cols, function(x) {return x.length;});
    var y_labels_max = d3.max(rows, function(x) {return x.length;});

    var XOffset = 10 + 7 * y_labels_max;  // +150;
    var YOffset = 10 + 6 * x_labels_max; // +100;

    if (!data || data.length == 0) {
        $(placeholder_id).html("Not enough data points to compare");
        return;
    }
    var X = cols.length;
    var Y = rows.length;

      // .offset(function(d) {
      //   if (col_accessions[d.x] == row_accessions[d.y]) {
      //     return [-30, 10]
      //   } else {
      //     return [-60, -10]
      //   }
      // })
    var tip = d3.tip()
      .attr('class', 'd3-tip')
      .offset([-9, 2])
      .html(function(d) {

          if (col_accessions[d.x] == row_accessions[d.y]) {
             var src = url_fingerprint.replace('0', col_accessions[d.x]);
             return 'Click to show <strong>' + rows[d.y] + '</strong><br><img height="45px" width="224px" style="padding: 10px;" src="' + src + '"</img>';
          }

          var src_x = url_fingerprint.replace('0', col_accessions[d.x]);
          var src_y = url_fingerprint.replace('0', row_accessions[d.y]);
          var value = +d.value || 0;
          var comparison = "Click to compare:<br><br> <strong>" + rows[d.y] + '</strong><br>';
          comparison += '<img height="45px" width="224px" height="45px" width="224px" style="padding: 10px;"  src="' + src_y + '"</img><br>';
          comparison += 'to <strong>' + cols[d.x] + "</strong><br>";
          comparison += '<img height="45px" width="224px" style="padding: 10px;"  src="' + src_x + '"</img><br><br>';
          // chisquare, hellinger, jensenshannon, cos, euclidean
          if (method == "chisquare") {
            return comparison + "Chi-square p-value = " + value;
          } else if (method == "cos") {
            return comparison + "Cosine distance = " + value;
          } else if (method == "hellinger") {
            return comparison + "Hellinger distance = " + value;
          } else if (method == "euclidean") {
            return comparison + "Euclidean distance = " + value;
          } else if (method == "jensenshannon") {
            return comparison + "Jensen-Shannon divergence = " + value;
          } else if (method == "chisqdist") {
            return comparison + "Chi-squared distance = " + value;
          } else if (method == "correlation") {
            return comparison + "Correlation coefficient = " + value;
          }
          return comparison;
      });

    var all_values = [];
    for(var i=0; i< data.length; i++)
    {
      if (data[i].x != data[i].y) {
         all_values.push(data[i].value);
      }
    }
    var maximum = Math.max.apply(Math, all_values);
    var minimum = Math.min.apply(Math, all_values);
    // colorScale = d3.scale.linear().domain([0.0, maximum]).range(["#00008B",  "#C80815"]);
    // colorScale = d3.scale.linear().domain([minimum, maximum]).range(["white",  "#C80815"]);
    // colorScale = d3.scale.linear().domain([minimum, 0.0, maximum]).range(["#FAFAFA", "#FAFAFA", "#904EAA"]);// "#9B58B5"]);
    // colorScale = d3.scale.linear().domain([minimum, 0.0, maximum]).range(["#FAFAFA", "#FAFAFA", "#9B58B5"]);  // gray - gray - violet


    // colorScale = d3.scale.linear().domain([minimum, 0.5*(maximum-minimum), maximum]).range(["#9B58B5", "yellow", "#ff8c00"]);
    // colorScale = d3.scale.linear().domain([0, 0.25*(maximum-minimum), 0.5*(maximum-minimum), 0.75*(maximum-minimum), maximum]).range(["#5e3c99","#b2abd2","#f7f7f7", "#fdb863", "#e66101"]);
    // colorScale = d3.scale.linear().domain([0, 0.9*(maximum-minimum), maximum]).range(["#b2abd2", "#fdb863", "#e66101"]);

    // colorScale = d3.scale.linear().domain([0, 0.9*(maximum-minimum), maximum]).range(["#ffffdf", "#fdae61", "#d7191c"]);    // yellow - orange - red
    // colorScale = d3.scale.linear().domain([0, minimum, 0.75*(maximum-0), maximum]).range(["#ffffdf", "#ffffdf", "#fdae61", "#d7191c"]);    // yellow - yellow - orange - red
    // colorScale = d3.scale.linear().domain([0, 0.9*(maximum-minimum), maximum]).range(["#ffffdf", "#b2abd2", "#904EAA"]); // yellow - violet - intence violet

    minimum = 0.0;
    // maximum = 1.0;

    var threshold = 0.5;
    // if (method == "hellinger") {
    //     threshold = 0.2;
    // } else if (method == "jensenshannon") {
    //     threshold = 0.08;
    // } else if (method == "cos") {
    //     threshold = 0.07;
    // }  else if (method == "chisquare") {
    //     threshold = 0.5;
    // } else if (method == "chisqdist") {
    //    threshold = 0.2;
    // }
    // threshold = Math.max(d3.quantile(all_values, 0.25),  d3.mean(all_values) - d3.deviation(all_values) * 3.0);
    // threshold = d3.quantile(all_values, 0.3);
    threshold =  d3.median(all_values) - d3.deviation(all_values) * 2.0;
    // colorScale = d3.scale.linear().domain([0, maximum]).range(["green", "#ffffef" ]);    // red - orange - white
    colorScale = d3.scale.linear().domain([0, threshold,  d3.mean(all_values), maximum]).range(["maroon", "red", "#ffffef" , "#ffffef" ]);    // red - orange - white
    // colorScale = d3.scale.linear().domain([0, maximum]).range(["#d7191c", "#ffffdf", ]);    // yellow - yellow - orange - red

    // if (method == "euclidean") {
    //     threshold = 0.001;
    //     colorScale = d3.scale.linear().domain([0, threshold, maximum]).range(["#d7191c",  "#fdae61", "#ffffdf",]);    // yellow - yellow - orange - red
    // }

       // "color": "#B3DA29"
       //  "color": "#78C0F2"
       //  "color": "#9B58B5"
       //  "color": "#E54C3C"
       //  "color": "#F2519D"

    $(placeholder_id).html("");
    var svg = d3.select(placeholder_id)
        .append("svg")
        .attr("width", w*X + spacer*(X-1) + XOffset + XOffset)
        .attr("height", h*Y + spacer*(Y-1) + YOffset);

    if (svg == null) {
      console.log("Could not initialize SVG in draw_matrix");
      return;
    }

    // init tooltip
    svg.call(tip);

    var plotCell = svg.selectAll(".map")
         // .selectAll(".cell")
         .data(data)
           // .enter().append("g")
           .enter().append("svg:rect")
             .attr('width', w)
             .attr('height',h)
             .attr('x', function(d) { 
                return d.x * (w + spacer) + XOffset;
             })
             .attr('y', function(d) {
                return d.y * (h + spacer) + YOffset;
             })
             .style('fill', function(d) {
                var value = +d.value || 0;
                return colorScale(value);
             })
             .on('click', function(d) {
                  if (col_accessions[d.x] == row_accessions[d.y]) {
                      var header = cols[d.x];
                      if (mutation_type_labels.columns.length > 0) {
                          header += ' (' + mutation_type_labels.columns + ')';
                      }
                      show_signature_modal(col_accessions[d.x], cols[d.x], header);
                  } else {
                      var header = "";
                      header += cols[d.x];
                      if (mutation_type_labels.columns.length > 0) {
                          header += ' (' + mutation_type_labels.columns + ')';
                      }
                      header += ' vs. '
                      header += rows[d.y];
                      if (mutation_type_labels.rows.length > 0) {
                          header += ' (' + mutation_type_labels.rows + ')';
                      }
                      show_comparison_modal(col_accessions[d.x], row_accessions[d.y], header);
                  }
             })
             .on('mouseover', function(d){
                  // tip.show(d, document.getElementById("heading1"));
                  // console.log(this);
                  d3.select(this).classed("cell-hover", true);
                  tip.show(d, this);
             })
             .on('mouseout', function(d) {
                  d3.select(".cell-hover").classed("cell-hover", false);
                  tip.hide(); 
              }); 

    //label rows
    var rowLabels = svg.append("g")
        .attr("class", "rowLabels")
        .selectAll(".rowLabel")
        .data(rows)
        .enter().append('text')
        .attr('y', function(d, i) {
            return YOffset + 2 + i*(h + spacer) + (h+spacer)/2;
        })
        .attr('x', XOffset - 5)
        .attr('class','matrix-label') 
        .attr('title', function(d) {return d;})
        .style('text-anchor','end')
        .text(function(d, i) {
            return d;
        });

    //label columns
    var colLabels = svg.append("g")
        .attr("class", "colLabels")
        .selectAll(".colLabel")
        .data(cols)
        .enter().append('text')
        .attr('x', 0)
        .attr('y', 0)
        .attr("transform", function(d, j) {
            return "translate("+ (XOffset + 2 + j*(w + spacer) + (w+spacer)/2) +", "+YOffset+") rotate(-45)"
        })
        .attr('class','matrix-label')
        .style('text-anchor','start')
        .text(function(d, i) {
            return d;
        });

}

function show_scatterplot(ref, uri_scatterplot, logscale){
      $.getJSON(uri_scatterplot, {logscale: logscale},  function(json) {
            draw_scatterplot(ref, json, logscale);

            $("#hist1_loading").text(json.xlabel);
            $("#hist2_loading").text(json.ylabel);
            $("#logplot_loading").text("");
            $("#logplot_label1").text(json.xlabel);
            $("#logplot_label2").text(json.ylabel);        
      });  
}


function show_similarity(ref, url){
      var $html = $(ref);
      $html.html("");
      $html.append($('<table>').addClass("table table-hover"));
      $html.append()
      $.getJSON(url, function(json) {
            // console.log(json);
            for (var i in json.data) {
                // two special cases:
                var d = json.data[i];
                // if (d.method == "chisquare") {
                //     var dt = $('<dt>').text(d.label);
                //     var dd = $('<dd>').html("<br>" +
                //       'Chi statistic: ' + d.value[0].toPrecision(3) + '<br>' + 
                //       'P-value: ' + d.value[1].toExponential(4));
                //     $(ref + " > dl").append(dt).append(dd);
                //     continue;
                // }
                // if (d.method == "euclidean") {
                //     var dt = $('<dt>').text(d.label);
                //     var dd = $('<dd>').html("<br>" + 
                //       'Distance: ' + d.value[0].toPrecision(3)); 
                //     $(ref + " > dl").append(dt).append(dd);
                //     continue;
                // }

                // the rest of the cases:

                var tr = $('<tr>');
                tr.append($('<td>').addClass("text-right").html(d.label));
                tr.append($('<td width="70%">').addClass("text-left").html(d.value[0].toPrecision(3)));
                // tr.append($('<td>').html(""));
                
                $(ref + " > table").append(tr);
            }
      });
}


function show_comparison_modal(id1, id2, header, data) {
      if (header) {
          $("#comparison_modal_header").text(header);
      } else {
          $("#comparison_modal_header").text("");
      }
      
      var uri_scatterplot = encodeURI(url_api_compare_pair_scatterplot.replace('s1', id1).replace('s2', id2));
      $("#scatterplot_logscale").off('change');
      $("#scatterplot_logscale").on('change', function(){
            show_scatterplot("#scatterplot", uri_scatterplot, this.checked);
      });

      d3.selectAll("#scatterplot" + " > *").remove();
      show_scatterplot(
          "#scatterplot",
          uri_scatterplot,
          $("#scatterplot_logscale").prop( "checked" ));

      d3.selectAll("#logplot" + " > *").remove();
      drawLogPlot(
          "#logplot",
          encodeURI(url_api_compare_pair_logplot.replace('s1', id1).replace('s2', id2)));

      d3.selectAll("#hist1" + " > *").remove();
      drawBarplot(
          "#hist1",
          encodeURI(url_signature_freqs.replace('0', id1)),
          max_width=750, max_height=280, click=undefined, show_legend=true);

      d3.selectAll("#hist2" + " > *").remove();
      drawBarplot(
          "#hist2",
          encodeURI(url_signature_freqs.replace('0', id2)),
          max_width=750, max_height=280, click=undefined, show_legend=false);
      
      show_similarity(
          "#similarity",
          encodeURI(url_api_compare_pair.replace('s1', id1).replace('s2', id2)));
      
      $("#comparisonModal").modal('show');
}

function scrollToAnchor(aid){
    var aTag = $("a[name='"+ aid +"']");
    $('html,body').animate({scrollTop: aTag.offset().top},'slow');
}

////////////////////////////////////////////////////////////////////////////////////
