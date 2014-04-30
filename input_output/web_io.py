"""
Andrew Hill
MacArthur Lab - 2014

Functions for handling input and output from the web.
"""

import urllib


def get_page_html(page_url):
    """
    Returns the html describing the page at the specified URL.

    @param page_url: URL to a specified website
    @type page_url: string
    @return: HTML describing the specified page
    @rtype: string
    @raise: IOError if an error code is returned
    @raise: ValueError if requested URL not found.
    """
    try:
        response_code = urllib.urlopen(page_url).code
    except IOError as e:
        raise ValueError('Requested URL not found. Please ensure http:// is included.')

    if response_code >= 400:
        raise IOError('Requested URL could not be reached. Response code: ' + response_code)

    return urllib.urlopen(page_url).read()