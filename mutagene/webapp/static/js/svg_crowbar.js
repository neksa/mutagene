// FROM http://bl.ocks.org/deanmalmgren/22d76b9c1f487ad1dde6
//
// apply the stylesheet to the svg to be sure to capture all of the stylings
//   applyStylesheets(svg_el)
//

// this is adapted (barely) from svg-crowbar
// https://github.com/NYTimes/svg-crowbar/blob/gh-pages/svg-crowbar-2.js#L211-L250
function applyStylesheets(svgEl) {

    // use an empty svg to compute the browser applied stylesheets
    var emptySvg = window.document.createElementNS("http://www.w3.org/2000/svg", 'svg');
    window.document.body.appendChild(emptySvg);
    var emptySvgDeclarationComputed = getComputedStyle(emptySvg);
    emptySvg.parentNode.removeChild(emptySvg);

    // traverse the element tree and explicitly set all stylesheet values
    // on an element. this is ripped from svg-crowbar
    var allElements = traverse(svgEl);
    var i = allElements.length;
    while (i--){
        explicitlySetStyle(allElements[i], emptySvgDeclarationComputed);
    }
}


function explicitlySetStyle (element, emptySvgDeclarationComputed) {
    var cSSStyleDeclarationComputed = getComputedStyle(element);
    var i, len, key, value;
    var computedStyleStr = "";
    for (i=0, len=cSSStyleDeclarationComputed.length; i<len; i++) {
        key=cSSStyleDeclarationComputed[i];
        value=cSSStyleDeclarationComputed.getPropertyValue(key);
        if (value!==emptySvgDeclarationComputed.getPropertyValue(key)) {
            computedStyleStr+=key+":"+value+";";
        }
    }
    element.setAttribute('style', computedStyleStr);
}


// traverse an svg and append all of the elements to the tree array. This
// ignores some elements that can appear in <svg> elements but whose
// children's styles should not be tweaked
function traverse(obj){
    var tree = [];
    var ignoreElements = {
        'script': undefined,
        'defs': undefined,
    };
    tree.push(obj);
    visit(obj);
    function visit(node) {
        if (node && node.hasChildNodes() && !(node.nodeName.toLowerCase() in ignoreElements)) {
            var child = node.firstChild;
            while (child) {
                if (child.nodeType === 1) {
                    tree.push(child);
                    visit(child);
                }
                child = child.nextSibling;
            }
        }
    }
    return tree;
}
