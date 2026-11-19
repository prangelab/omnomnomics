(function () {
  var colorBinding = new Shiny.InputBinding();

  function showValue(element) {
    element.parentElement.querySelector('.explorer-color-hex').textContent = element.value;
  }

  $.extend(colorBinding, {
    find: function (scope) {
      return $(scope).find('input.explorer-color-input');
    },
    getValue: function (element) {
      return element.value;
    },
    subscribe: function (element, callback) {
      $(element).on('change.explorerColors', function () {
        showValue(element);
        callback(false);
      });
    },
    unsubscribe: function (element) {
      $(element).off('.explorerColors');
    },
    receiveMessage: function (element, message) {
      if (message.value !== undefined) {
        element.value = message.value;
        showValue(element);
        $(element).trigger('change');
      }
    }
  });

  Shiny.inputBindings.register(colorBinding, 'omnomnomics.groupColor');
})();
