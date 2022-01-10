export class ValueSelector extends HTMLElement {
    constructor() {
        super();

        this.multiple = this.hasAttribute("multiple");
        let defaultHeader = this.hasAttribute("header")? this.getAttribute("header"): "Chosen genes";
        let root = this.attachShadow({ mode: 'open'});

        let style = document.createElement("style");
        style.textContent = `

`;
        root.appendChild(style);

        // Apply external styles to the shadow dom
        const cssShard = document.createElement('link');
        cssShard.setAttribute('rel', 'stylesheet');
        cssShard.setAttribute('href', '/styles/shards-dashboards.1.1.0.min.css');
        root.appendChild(cssShard);
        // Apply external styles to the shadow dom
        const cssApp = document.createElement('link');
        cssApp.setAttribute('rel', 'stylesheet');
        cssApp.setAttribute('href', '/styles/app.css');
        root.appendChild(cssApp);

        if (this.hasAttribute("stylesheet")) {
            let extcss = document.createElement("link");
            extcss.setAttribute("rel", "stylesheet");
            extcss.setAttribute("href", this.getAttribute("stylesheet"));
            root.appendChild(extcss);
        }
        
        let text = document.createElement('input');
        text.setAttribute('type', 'text');
        if (this.hasAttribute("placeholder")) {
            text.setAttribute("placeholder", this.getAttribute("placeholder"));
        }

        text.setAttribute("class", "form-control");

        let filters = document.createElement("fieldset");

        let endpointSelector = document.createElement("select");
        endpointSelector.setAttribute("class", "form-control");

        for (let e of this.getElementsByTagName("endpoint")) {
            let o = document.createElement("option");
            o.setAttribute("value", e.getAttribute("name"));
            o.textContent = e.hasAttribute("display-text")? e.getAttribute("display-text"): e.textContent;
            endpointSelector.appendChild(o);
        }

        endpointSelector.firstChild.setAttribute("selected","");

        endpointSelector.addEventListener("change", c => {
            this.find(".holder").classList.add("hidden");
            this.find("ul.results").innerHTML = "";
            this.find("fieldset").innerHTML = ""; // clear old filters
            let endpoint = this.querySelector(`endpoint[name="${c.currentTarget.value}"]`);
            for (let f of endpoint.querySelectorAll("filter")) {
                if (f.hasAttribute("label")) {
                    let l = document.createElement("label");
                    l.textContent = f.getAttribute("label");
                    filters.appendChild(l);
                }
                if (f.getAttribute("type") == "list") {
                    let s = document.createElement("select");
                    s.setAttribute("name", f.getAttribute("parameter"));
                    for (let e of f.querySelectorAll("entry")) {
                        let o = document.createElement("option");
                        o.setAttribute("value", e.getAttribute("value"));
                        o.textContent = e.hasAttribute("display-text")? e.getAttribute("display-text"): e.getAttribute("value");
                        if (e.hasAttribute("default")) o.setAttribute("selected", "");
			s.appendChild(o);
                    }
                    filters.appendChild(s);
                }
            }
        })

        if(this.getElementsByTagName("endpoint").length == 1) {
            endpointSelector.classList.add("hidden");
        }

        let buttons = "";
        if (this.hasAttribute("buttons")) {
            buttons = document.createElement("div");
            buttons.setAttribute("id", "buttons");
            
            let buttonlist = this.getAttribute("buttons").split(":");

            if (buttonlist.includes("remove") || buttonlist.includes("all")) {
                let removeButton = document.createElement("button");
                removeButton.textContent = "Remove All";
                removeButton.addEventListener("click", c => {
                    c.stopPropagation();
                    this.removeAll();
                    this.find(".header").textContent = defaultHeader;
                    this.dispatchEvent(new Event("update"));
                });
                buttons.appendChild(removeButton);
            }
            if (buttonlist.includes("select") || buttonlist.includes("all")) {
                let button = document.createElement("button");
                button.textContent = "Select All";
                button.addEventListener("click", c => {
                    c.stopPropagation();
                    this.findAll("ul.choices li:not(.selected)").forEach(i => i.classList.add("selected"));
                    this.dispatchEvent(new Event("update"));
                });
                buttons.appendChild(button);
            }
            if (buttonlist.includes("deselect") || buttonlist.includes("all")) {
                let button = document.createElement("button");
                button.textContent = "Deselect All";
                button.addEventListener("click", c => {
                    c.stopPropagation();
                    this.findAll("ul.choices li.selected").forEach(i => i.classList.remove("selected"));
                    this.dispatchEvent(new Event("update"));
                });
                buttons.appendChild(button);
            }
        }
            
        let holder = document.createElement("div");
        holder.classList.add("holder");
        holder.classList.add("hidden");

        let msg = document.createElement("label");

        let close = document.createElement("label");
        close.textContent = "Close list";
        close.classList.add("close-box");

        close.addEventListener("click", c => {
            c.stopPropagation();
            holder.classList.add("hidden");
        })

        let list = document.createElement("ul");
        list.classList.add("results");
        
        list.addEventListener("click", c => {
            c.stopPropagation();
            this.setHeader(defaultHeader);

            let src = c.target;
            while (src.tagName.toLowerCase() != "li") src = src.parentNode;

            if (src.dataset.action == "replace") {
                this.removeAll();
                holder.classList.add("hidden");
                this.setHeader(src.textContent);
            }

            if (src.dataset.list) {
                this.push(src.dataset.list.split(","));
            } else {
                this.push(src.textContent);
            }

            if (src.dataset.action != "replace") src.remove();

            if (list.childNodes.length == 0) holder.classList.add("hidden");
        })

        holder.append(close,document.createElement("br"),msg,list);

        let header = document.createElement("h4");
        header.classList.add("header");
        header.textContent = defaultHeader;

        let choices = document.createElement("ul");
        choices.classList.add("choices");

        choices.addEventListener("click", c => {
            c.stopPropagation();

            let src = c.target;

            while (src.tagName.toLowerCase() != "li") src = src.parentNode;

            if (this.multiple) {
                if (src.classList.contains("selected")) {
                    src.classList.remove("selected");
                } else {
                    src.classList.add("selected");
                }
            } else {
                if (src.classList.contains("selected")) {
                    // allow deselection of selected item in single mode;
                    // will trigger update event with no selected items.
                    src.classList.remove("selected");
                } else {
                    if (choices.querySelector("li.selected")) choices.querySelector("li.selected").classList.remove("selected");
                    src.classList.add("selected");
                }
            }
    
            this.dispatchEvent(new Event("update"));
    
        })

        let controller = new AbortController();

        text.addEventListener('input', ti => {
            let endpoint = this.querySelector(`endpoint[name="${endpointSelector.value}"]`);
            
            // if this endpoint is triggered by something other than text input, check we've received a 'proper' InputEvent** & bail out.
            // **Endpoints with non-input triggers operate by dispatching an Event("input") on text, which will not have an inputType property.
            if (ti.inputType && endpoint.hasAttribute("trigger") && endpoint.getAttribute("trigger") != "input") return;

            controller.abort(); // cancel any ongoing fetch from this handler
    
            // set up the search url
            let searchterm = ti.target.value;

            if (searchterm.length == 0 && !endpoint.hasAttribute("searchblank")) {
                // do not search for an empty string! clear and hide list
                while(list.firstChild) list.removeChild(list.lastChild);
                holder.classList.add("hidden");
                return;
            }

            let additionalterms = {};
            if (searchterm.length > 0) {
                for (let f of filters.querySelectorAll("select, input")) {
                    additionalterms[f.getAttribute("name")] = f.value;
                }
            }

            let getData;

            if (endpoint.hasAttribute("url")) {
                let connector = endpoint.getAttribute("url").includes('?')? "&": "?";

                let url = endpoint.getAttribute("url").concat(connector, endpoint.getAttribute("parameter"), '=', searchterm);
                
                Object.keys(additionalterms).forEach(k => url = url.concat("&", k, "=", additionalterms[k]))
                // reset the controller
                controller = new AbortController();
        
                getData = fetch(url, { signal: controller.signal })
                    .then(response => response.json())

            } else if (endpoint.hasAttribute("function")) {
                let datafunction = window[endpoint.getAttribute("function")];
                if (!datafunction || typeof(datafunction) != "function") throw new ReferenceError(`Data source "${endpoint.getAttribute("function")}" is not a known function in the global scope.`);

                getData = new Promise((resolve, reject) => {
                    try {
                        resolve(datafunction(searchterm, additionalterms));
                    } catch(e) {
                        reject(e);
                    }
                })
            } else {
                throw `No source-type attribute (url|function) in endpoint "${endpoint.getAttribute("name")}"}`
            }
                
            getData.then(data => {
                    while(list.firstChild) list.removeChild(list.lastChild);

                    if (endpoint.hasAttribute("preprocessor")) {
                        let preprocessor = window[endpoint.getAttribute("preprocessor")];
                        if (typeof(preprocessor) == "function") data = preprocessor(data);
                        else throw new ReferenceError(`Data preprocessor "${endpoint.getAttribute("preprocessor")}" is not a known function in the global scope.`)
                    } 

                    let results = data[endpoint.getAttribute("property")];

                    if (results.length == 0 && endpoint.getAttribute("empty") == "noaction") return;

                    if (endpoint.getAttribute("mode") == 'list') {
                        results = results.filter(item => !this.items.includes(item)).sort();
                        

                        results.forEach(i => { 
                            let e = document.createElement("li");
                            e.textContent = i;
                            list.appendChild(e);
                        });
                    } else if (endpoint.getAttribute("mode") == "groups") {
                        Object.keys(results).forEach(i => {
                            let e = document.createElement("li");
                            e.textContent = i;
                            e.setAttribute("data-list", results[i].join(","));
                            e.setAttribute("data-count", results[i].length);
                            list.appendChild(e);
                        })
                    }

                    if(list.childNodes.length == 0) {
                        let emptymessage = "No matching items";
                        if (endpoint.hasAttribute("empty") && endpoint.getAttribute("empty").startsWith("msg"))
                            emptymessage = endpoint.getAttribute("empty").split(":")[1];

                        msg.textContent = emptymessage;
                    } else {
                        msg.textContent = "";
                    }

                    if (endpoint.getAttribute("action") == "replace") {
                        list.childNodes.forEach(n => n.setAttribute("data-action", "replace"));
                    }
			
		            holder.classList.remove("hidden");
                })
                .catch(e => {
                    if (e.name == "AbortError") {
                        // pass - we have cancelled the fetch by re-entering the handler
                    } else {
                        console.log(e);
                    }
                });
        });
        text.addEventListener("focus", fe => {
            if(list.childNodes.length > 0) {
                if (list.firstChild.textContent.toLowerCase().includes(text.value.toLowerCase())) {
                    holder.classList.remove("hidden");
                } else {
                    while(list.firstChild) list.removeChild(list.lastChild);
                }
            } else if(text.value.length == 0 && this.querySelector(`endpoint[name="${endpointSelector.value}"]`).hasAttribute("searchblank")) {
                text.dispatchEvent(new Event("input"));
            } 
        })
        text.addEventListener("keydown", ke => {
            let endpoint = this.querySelector(`endpoint[name="${endpointSelector.value}"]`);

            // handle case of endpoint that is triggered by a key event (usually enter):
            if (endpoint.hasAttribute("trigger") && endpoint.getAttribute("trigger").toLowerCase() == ke.key.toLowerCase()) {
                text.dispatchEvent(new Event("input"));
                return;
            }

            // by default handle a keypress of 'enter' to simulate a click if any entry in 'results' exactly matches the current value of text
            if (ke.key.toLowerCase() == "enter") {
                Array.from(this.findAll("ul.results li")).filter(e => e.textContent == ke.currentTarget.value).forEach(e => e.dispatchEvent(new MouseEvent('click', {
                    view: window,
                    bubbles: true,
                    cancelable: true
                  })));
                return;
            }
        })
        root.append(endpointSelector,text,filters,buttons,holder,header,choices);

        if(this.dataset.init) {
            this.push(this.dataset.init.split(",").map(s => s.trim()))
        }
    }

    get items() {
        return Array.from(this.shadowRoot.querySelectorAll("ul.choices li label"), l => l.textContent);
    }

    get selectedItems() {
        return Array.from(this.shadowRoot.querySelectorAll("ul.choices li.selected label"), l => l.textContent);
    }

    push(value) {
        if (Array.isArray(value)) {
            // add multiple items to the selection
            value.forEach(val => this.add(val));

        } else {
            // add one item to the selection
            this.add(value);
        }

        this.dispatchEvent(new Event("update"));
    }

    add(value) {
        if (!value) {
            console.log(`Attempt to add empty, null or undefined value.`);
            return;
        } else if (this.items.includes(value)) {
            console.log(`Attempt to re-add existing item "${ value }".`);
            return;
        }
        
        let item = document.createElement("li");
        const choices = this.find("ul.choices");
        
        let label = document.createElement("label");
        label.textContent = value;

        let close = document.createElement("div");
        close.setAttribute("class", "close");
        close.setAttribute("data-value", value);
        close.textContent = "x";

        close.addEventListener("click", c => { 
            c.stopPropagation();

            const src = c.target;
            const data = src.dataset.value;
            const searchterm = this.find("input[type='text']").value;
            
            if (searchterm.length > 0 && data.includes(searchterm)) {
                const list = this.find("ul.results");
                let li = document.createElement("li");
                li.textContent = data;
                let i=0;
                while(list.childNodes[i] && list.childNodes[i].textContent < data) i++;
                list.insertBefore(li, list.childNodes[i]);
            }
            src.parentNode.remove();

            this.dispatchEvent(new Event("update"));
        });

        item.append(label, close);
        if (this.multiple) item.classList.add("selected");
        choices.appendChild(item);
        return item;
    }

    select(item) {
        this.find(`button[data-value="${item}"]`).parentNode.classList.add("selected");
    }

    deselect(item) {
        this.find(`button[data-value="${item}"]`).parentNode.classList.remove("selected");
    }

    remove(item) {
        this.find(`.choices li:has(button[data-value='${item}'])`).remove();
    }

    removeAll() {
        this.find(".choices").innerHTML = "";
    }

    setHeader(text) {
        this.find(".header").textContent = text;
    }

    find(selector) {
        return this.shadowRoot.querySelector(selector);
    }

    findAll(selector) {
        return this.shadowRoot.querySelectorAll(selector);
    }
}
