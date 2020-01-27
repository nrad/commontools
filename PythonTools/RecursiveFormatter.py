"""
https://gist.github.com/26f09ae6f551d12e4249.git
"""

from string import Formatter


class RecursiveFormatter(Formatter):
  def _contains_underformatted_field_names(self, format_tuple):
    literal_text, field_name, format_spec, conversion = format_tuple
    if field_name is not None:
      return any(component[1] is not None for component in self.parse(field_name))
    return False

  def _split_format_tuple(self, format_tuple):
    literal_text, field_name, format_spec, conversion = format_tuple
    
    yield (literal_text + '{', None, None, None)
    
    for format_tuple in self.parse(field_name):
      yield format_tuple
    
    literal_text = ''.join([
        '' if not conversion else ('!' + conversion),
        '' if not format_spec else (':' + format_spec),
        '}'])

    yield (literal_text, None, None, None)

  def parse(self, format_string):
    def iter_tuples():
      for format_tuple in super(RecursiveFormatter, self).parse(format_string):
        if not self._contains_underformatted_field_names(format_tuple):
          yield format_tuple
        else:
          for subtuple in self._split_format_tuple(format_tuple):
            yield subtuple
    return list(iter_tuples())

  def vformat(self, format_string, args, kwargs):
    while True:
      format_string = super(RecursiveFormatter, self).vformat(format_string, args, kwargs)

      # if any field_names remain uninterpolated
      if any(format_tuple[1] is not None for format_tuple in self.parse(format_string)):
        continue
      else:
        return format_string


